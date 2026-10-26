"""Numba engine for the per-node merge rules (R1-R3, RD and quartet tests).

Same algorithm as the pure-Python reference in :func:`ancestra.merge.decide_node`
(``engine="python"``); data structures are CSR arrays (neighbours sorted by unit id,
binary-search lookups) and an array union-find with per-component member lists.
Tie-breaking can differ from the Python reference only between exactly equal scores
(Python uses dict insertion order, numba uses the smaller unit / gene id);
tests/test_merge.py checks agreement on the toy simulations.
"""
from __future__ import annotations

import numba as nb
import numpy as np

NEG = -np.inf


@nb.njit(cache=True)
def csr_undirected(n, a, b, s):
    deg = np.zeros(n + 1, np.int64)
    for i in range(len(a)):
        deg[a[i] + 1] += 1
        deg[b[i] + 1] += 1
    ptr = np.cumsum(deg)
    nbr = np.empty(ptr[-1], np.int64)
    sc = np.empty(ptr[-1], np.float64)
    fill = ptr[:-1].copy()
    for i in range(len(a)):
        x = a[i]
        y = b[i]
        nbr[fill[x]] = y
        sc[fill[x]] = s[i]
        fill[x] += 1
        nbr[fill[y]] = x
        sc[fill[y]] = s[i]
        fill[y] += 1
    for u in range(n):
        lo = ptr[u]
        hi = ptr[u + 1]
        if hi - lo > 1:
            o = np.argsort(nbr[lo:hi], kind="mergesort")
            nbr[lo:hi] = nbr[lo:hi][o]
            sc[lo:hi] = sc[lo:hi][o]
    return ptr, nbr, sc


@nb.njit(cache=True)
def csr_directed_sorted(n, a, b, s):
    """a sorted ascending, and b sorted within equal a (as produced by _aggregate)."""
    ptr = np.zeros(n + 1, np.int64)
    for i in range(len(a)):
        ptr[a[i] + 1] += 1
    ptr = np.cumsum(ptr)
    return ptr, b.astype(np.int64), s.astype(np.float64)


@nb.njit(cache=True)
def lookup(ptr, nbr, sc, u, v):
    lo = ptr[u]
    hi = ptr[u + 1]
    if hi <= lo:
        return NEG
    i = lo + np.searchsorted(nbr[lo:hi], v)
    if i < hi and nbr[i] == v:
        return sc[i]
    return NEG


@nb.njit(cache=True)
def find(parent, x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


@nb.njit(cache=True)
def union(parent, size, head, tail, nxt, a, b):
    ra = find(parent, a)
    rb = find(parent, b)
    if ra == rb:
        return ra
    if size[ra] < size[rb]:
        ra, rb = rb, ra
    parent[rb] = ra
    size[ra] += size[rb]
    nxt[tail[ra]] = head[rb]
    tail[ra] = tail[rb]
    return ra


@nb.njit(cache=True)
def _dlog(v):
    if v < 1e-6:
        v = 1e-6
    if v > 1.0:
        v = 1.0
    return -np.log(v)


@nb.njit(cache=True)
def rel_div(optr, onbr, osc, a, b, s_ab):
    la = optr[a + 1] - optr[a]
    lb = optr[b + 1] - optr[b]
    if la == 0 or lb == 0:
        return np.nan
    x, y = a, b
    if la > lb:
        x, y = b, a
    best_v = -1.0
    va_b = 0.0
    vb_b = 0.0
    found = False
    for i in range(optr[x], optr[x + 1]):
        z = onbr[i]
        vy = lookup(optr, onbr, osc, y, z)
        if vy != NEG:
            if osc[i] + vy > best_v:
                best_v = osc[i] + vy
                va_b = osc[i]
                vb_b = vy
                found = True
    if not found:
        return np.nan
    den = 0.5 * (_dlog(va_b) + _dlog(vb_b))
    if den <= 0:
        return np.nan
    return _dlog(s_ab) / den


@nb.njit(cache=True)
def rd_values(optr, onbr, osc, ca, cb, cs, idx):
    out = np.full(len(idx), np.nan)
    for j in range(len(idx)):
        e = idx[j]
        out[j] = rel_div(optr, onbr, osc, ca[e], cb[e], cs[e])
    return out


@nb.njit(cache=True)
def _has_child_count(head, nxt, child, root, c, mark, stamp):
    """Return (has member of child c, number of distinct children) of component root."""
    has = False
    nd = 0
    p = head[root]
    while p != -1:
        if mark[child[p]] != stamp:
            mark[child[p]] = stamp
            nd += 1
        if child[p] == c:
            has = True
        p = nxt[p]
    return has, nd


@nb.njit(cache=True)
def rules(n, k, child, ca, cb, cs, sa, sb, ss, idx, veto, min_score,
          inpar_on, inpar_factor, rule_code, rd_inpar, rd_factor, med_r,
          repair_on, repair_ratio, optr, onbr, osc, has_out, cs_r3, sup, veto_on, veto_min):
    """R1 (greedy mutual-best union with in-paralog constraint), R2 (absorption),
    R3 (fission repair). rule_code: 0 ortholog, 1 self, 2 either, 3 quartet.
    Returns labels and counters [r1_unions, r1_blocked, r2_absorbed, r2_rd_absorbed,
    r2_quartet_tests, r3_repaired]."""
    cnt = np.zeros(6, np.int64)
    keep = cs > min_score
    xptr, xnbr, xsc = csr_undirected(n, ca[keep], cb[keep], cs[keep])
    sptr, snbr, ssc = csr_undirected(n, sa, sb, ss)
    # context evidence (synteny): R3 selection scores and support, same CSR layout as xsc
    _p2, _n2, xsel = csr_undirected(n, ca[keep], cb[keep], cs_r3[keep])
    _p3, _n3, xsup = csr_undirected(n, ca[keep], cb[keep], sup[keep])
    maxsup = np.zeros(n)
    for u in range(n):
        for i in range(xptr[u], xptr[u + 1]):
            if xsup[i] > maxsup[u]:
                maxsup[u] = xsup[i]
    best_any = np.full(n, NEG)
    for i in range(len(ca)):
        if cs[i] > best_any[ca[i]]:
            best_any[ca[i]] = cs[i]
        if cs[i] > best_any[cb[i]]:
            best_any[cb[i]] = cs[i]
    parent = np.arange(n)
    size = np.ones(n, np.int64)
    head = np.arange(n)
    tail = np.arange(n)
    nxt = np.full(n, -1, np.int64)
    markP = np.full(k, -1, np.int64)
    maxsim = np.full(k, NEG)
    mark = np.full(k, -1, np.int64)
    stamp = 0
    # ---------------- R1
    for j in range(len(idx)):
        e = idx[j]
        if veto[j]:
            continue
        x = ca[e]
        y = cb[e]
        s = cs[e]
        rx = find(parent, x)
        ry = find(parent, y)
        if rx == ry:
            continue
        stamp += 1
        p = head[rx]
        while p != -1:
            markP[child[p]] = stamp
            maxsim[child[p]] = NEG
            p = nxt[p]
        common = False
        q = head[ry]
        while q != -1:
            c = child[q]
            if markP[c] == stamp:
                common = True
                p = head[rx]
                while p != -1:
                    if child[p] == c:
                        v = lookup(sptr, snbr, ssc, p, q)
                        if v > maxsim[c]:
                            maxsim[c] = v
                    p = nxt[p]
            q = nxt[q]
        ok = True
        if common:
            q = head[ry]
            while q != -1:
                c = child[q]
                if markP[c] == stamp and maxsim[c] < inpar_factor * s:
                    ok = False
                    break
                q = nxt[q]
        if ok:
            union(parent, size, head, tail, nxt, x, y)
            cnt[0] += 1
        else:
            cnt[1] += 1
    # ---------------- R2
    if inpar_on:
        for _pass in range(5):
            changed = 0
            cu_l = []
            cs_l = []
            cv_l = []
            for u in range(n - 1, -1, -1):       # descending u: ties in s keep u-descending order
                if xptr[u + 1] > xptr[u] and size[find(parent, u)] == 1:
                    bi = xptr[u]
                    for i in range(xptr[u], xptr[u + 1]):
                        if xsc[i] > xsc[bi]:
                            bi = i
                    cu_l.append(u)
                    cs_l.append(xsc[bi])
                    cv_l.append(xnbr[bi])
            m = len(cu_l)
            if m == 0:
                break
            ua = np.array(cu_l)
            sa_ = np.array(cs_l)
            va = np.array(cv_l)
            order = np.argsort(-sa_, kind="mergesort")
            for jj in range(m):
                t = order[jj]
                u = ua[t]
                s = sa_[t]
                v = va[t]
                if size[find(parent, u)] != 1:
                    continue
                g = find(parent, v)
                cu = child[u]
                stamp += 1
                has, nd = _has_child_count(head, nxt, child, g, cu, mark, stamp)
                if not has or nd < 2:
                    continue
                ok = False
                w = head[g]
                while w != -1:
                    if child[w] != cu:
                        w = nxt[w]
                        continue
                    sw = lookup(sptr, snbr, ssc, u, w)
                    if sw == NEG:
                        w = nxt[w]
                        continue
                    if rule_code == 3 and has_out:
                        swv = lookup(xptr, xnbr, xsc, w, v)
                        bo = -1.0
                        ou = 0.0
                        ov = 0.0
                        fo = False
                        if swv != NEG:
                            for i in range(optr[w], optr[w + 1]):
                                z = onbr[i]
                                a_u = lookup(optr, onbr, osc, u, z)
                                a_v = lookup(optr, onbr, osc, v, z)
                                if a_u != NEG and a_v != NEG and osc[i] + a_v > bo:
                                    bo = osc[i] + a_v
                                    ou = a_u
                                    ov = a_v
                                    fo = True
                        if fo:
                            cnt[4] += 1
                            if sw * ov >= inpar_factor * swv * ou:
                                ok = True
                                break
                            w = nxt[w]
                            continue
                    if (rule_code == 1 or rule_code == 2) and sw >= inpar_factor * s:
                        ok = True
                        break
                    if rule_code == 0 or rule_code == 2 or rule_code == 3:
                        dbest = NEG
                        for i in range(xptr[w], xptr[w + 1]):
                            if find(parent, xnbr[i]) == g and xsc[i] > dbest:
                                dbest = xsc[i]
                        if sw >= inpar_factor * dbest:
                            ok = True
                            break
                    if rd_inpar and not np.isnan(med_r):
                        r = rel_div(optr, onbr, osc, u, w, sw)
                        if not np.isnan(r) and r <= rd_factor * med_r:
                            ok = True
                            cnt[3] += 1
                            break
                    w = nxt[w]
                if ok and veto_on and maxsup[u] >= veto_min:
                    if lookup(_p3, _n3, xsup, u, v) <= 0:
                        ok = False          # synteny precision veto (support 0 vs >= veto_min)
                if ok:
                    union(parent, size, head, tail, nxt, u, v)
                    cnt[2] += 1
                    changed += 1
            if changed == 0:
                break
    # ---------------- R3
    if repair_on:
        pg_val = np.full(n, NEG)
        hc_mark = np.full(n, -1, np.int64)
        hc_val = np.zeros(n, np.bool_)
        ub_g = np.full(n, -1, np.int64)
        ub_s = np.full(n, NEG)
        ustamp = 0
        for _it in range(20):
            singles = []
            for u in range(n):
                if xptr[u + 1] > xptr[u] and size[find(parent, u)] == 1:
                    singles.append(u)
            if len(singles) == 0:
                break
            qb_u = nb.typed.Dict.empty(nb.types.int64, nb.types.int64)
            qb_s = nb.typed.Dict.empty(nb.types.int64, nb.types.float64)
            ub_list = []
            for u in singles:
                ub_g[u] = -1
                cu = child[u]
                ustamp += 1
                touched = []
                for i in range(xptr[u], xptr[u + 1]):
                    g = find(parent, xnbr[i])
                    if hc_mark[g] != ustamp:
                        hc_mark[g] = ustamp
                        stamp += 1
                        has, nd = _has_child_count(head, nxt, child, g, cu, mark, stamp)
                        hc_val[g] = has
                        pg_val[g] = NEG
                        if not has:
                            touched.append(g)
                    if hc_val[g]:
                        continue
                    if xsel[i] > pg_val[g]:
                        pg_val[g] = xsel[i]
                if len(touched) == 0:
                    continue
                bg = touched[0]
                for g in touched:
                    if pg_val[g] > pg_val[bg]:
                        bg = g
                s = pg_val[bg]
                raw = NEG
                for i in range(xptr[u], xptr[u + 1]):
                    if find(parent, xnbr[i]) == bg and xsc[i] > raw:
                        raw = xsc[i]
                if raw < repair_ratio * best_any[u]:
                    continue
                ub_g[u] = bg
                ub_s[u] = s
                ub_list.append(u)
                for g in touched:
                    key = g * k + cu
                    if key not in qb_s or pg_val[g] > qb_s[key]:
                        qb_s[key] = pg_val[g]
                        qb_u[key] = u
            merged = 0
            if len(ub_list) > 0:
                ul = np.array(ub_list)
                sv = np.empty(len(ul))
                for i in range(len(ul)):
                    sv[i] = -ub_s[ul[i]]
                order = np.argsort(sv, kind="mergesort")
                for jj in range(len(ul)):
                    u = ul[order[jj]]
                    g = ub_g[u]
                    s = ub_s[u]
                    if size[find(parent, u)] != 1 or find(parent, g) != g:
                        continue
                    key = g * k + child[u]
                    if key not in qb_u or qb_u[key] != u:
                        continue
                    stamp += 1
                    has, nd = _has_child_count(head, nxt, child, g, child[u], mark, stamp)
                    if has:
                        continue
                    if size[g] == 1:
                        q = head[g]
                        if ub_g[q] != -1 and ub_g[q] != find(parent, u):
                            continue
                        if lookup(xptr, xnbr, xsc, u, q) < repair_ratio * best_any[q]:
                            continue
                    if veto_on and maxsup[u] >= veto_min:
                        bs = 0.0
                        for i in range(xptr[u], xptr[u + 1]):
                            if find(parent, xnbr[i]) == g and xsup[i] > bs:
                                bs = xsup[i]
                        if bs <= 0:
                            continue
                    union(parent, size, head, tail, nxt, u, g)
                    cnt[5] += 1
                    merged += 1
            for u in ub_list:
                ub_g[u] = -1
            if merged == 0:
                break
    labels = np.empty(n, np.int64)
    for i in range(n):
        labels[i] = find(parent, i)
    return labels, cnt
