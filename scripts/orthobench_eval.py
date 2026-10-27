#!/usr/bin/env python3
"""OrthoBench evaluation helper: official scores + per-RefOG TP/FP/FN + dev/test subsets + paired bootstrap.

Scoring is the official Open_Orthobench `BENCHMARKS/benchmark.py` (unmodified, imported):
  * the readers (`read_orthogroups_smart`, `read_orthogroups`), `check_orthogroups`,
    `read_refogs`, `read_uncertain_refogs` are the upstream functions;
  * the headline F/precision/recall/#exact are obtained by calling upstream
    `calculate_benchmarks_pairwise` (reported as "official");
  * the per-RefOG counts are a re-implementation of the body of that function's loop (same
    uncertain-gene removal, same n-choose-2 pair counting, same 1/(n-1) RefOG normalisation, q_even=True),
    returning thisTP/N, thisFP/N, thisFN/N for every RefOG. Summing over all 70 RefOGs in the same
    order as upstream reproduces the official F/precision/recall to floating-point identity; this is
    asserted at every run (`verify` block in the JSON, `--no-verify` to skip upstream call).
  * dev/test subsets (splits/orthobench_split.json, root repo) = F/P/R from the SUM of the per-RefOG
    normalised TP/FP/FN over the subset (this is exactly the upstream formula restricted to the subset).

CLI (backward compatible with run_orthobench.sh):
  orthobench_eval.py <BENCHMARKS dir> <orthogroups file> [--basic] [--json out.json] [--label NAME]
                     [--split splits/orthobench_split.json] [--per-refog-tsv out.tsv]
  orthobench_eval.py --compare A.json B.json [--subset all|dev|test] [-B 2000] [--seed 11] [--ni-margin 3]
     paired bootstrap over RefOGs of dF = F(A) - F(B) (JSONs written by the first form). With --ni-margin M
     (percentage points) a non-inferiority verdict is added: "comparable" iff the lower 95% bound >= -M
     (DESIGN §5 of the root repo: M = 3 pp = 0.03 F as a fraction).
Python API: score_orthogroups(...), subset_scores(per_refog, ids), bootstrap_delta_f(perA, perB, ids, B=2000, seed=11).
Units: F/precision/recall in percent (as upstream); `F_frac` etc. are fractions.
"""
import argparse, contextlib, importlib.util, io, json, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SPLIT = os.path.abspath(os.path.join(HERE, "..", "..", "splits", "orthobench_split.json"))


def load_official(bench_dir):
    spec = importlib.util.spec_from_file_location("ob", os.path.join(bench_dir, "benchmark.py"))
    ob = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ob)
    ob.__file__ = os.path.join(bench_dir, "benchmark.py")  # get_expected_genes uses dirname(__file__)
    return ob


def per_refog_counts(ref_ogs, uncert_genes, pred_ogs):
    """Per-RefOG normalised (TP, FP, FN) exactly as in upstream calculate_benchmarks_pairwise(q_even=True,
    q_remove_uncertain=True). Returns list of dicts (index i -> RefOG%03d, i+1)."""
    gene2og = {}
    for j, og in enumerate(pred_ogs):
        for g in og:
            gene2og.setdefault(g, []).append(j)
    out = []
    for i, (refOg, uncert) in enumerate(zip(ref_ogs, uncert_genes)):
        # genes of the full RefOG determine candidate OGs (upstream tests overlap with the full RefOG first)
        cand = sorted({j for g in refOg for j in gene2og.get(g, ())})
        refOg = refOg.difference(uncert)
        nRefOG = len(refOg)
        not_present = set(refOg)
        tp = fp = fn = 0.0
        nsplit = 0
        for j in cand:
            predOg = pred_ogs[j].difference(uncert)
            overlap = len(refOg.intersection(predOg))
            if overlap > 0:
                nsplit += 1
                not_present = not_present.difference(predOg)
                tp += overlap * (overlap - 1) / 2
                fp += overlap * (len(predOg) - overlap)
                fn += (nRefOG - overlap) * overlap
        fn += len(not_present) * (nRefOG - 1)
        assert fn % 2 == 0
        fn /= 2
        N = float(nRefOG - 1)
        out.append(dict(refog="RefOG%03d" % (i + 1), n_genes=nRefOG, n_pairs=nRefOG * (nRefOG - 1) // 2,
                        TP_raw=tp, FP_raw=fp, FN_raw=fn, n_pred_fragments=nsplit,
                        exact=bool(fn == 0 and fp == 0),
                        TP=tp / N, FP=fp / N, FN=fn / N))
    return out


def prf(tp, fp, fn):
    """(F, precision, recall) in percent; zeros when undefined (TP=0), matching run_orthobench convention."""
    if tp <= 0:
        return 0.0, 0.0, 0.0
    p = tp / (tp + fp)
    r = tp / (tp + fn)
    f = 2 * p * r / (p + r)          # same operation order as upstream (bit-identical floats)
    return 100. * f, 100. * p, 100. * r


def subset_scores(per_refog, ids=None):
    """Sum per-RefOG normalised TP/FP/FN over `ids` (iterable of 'RefOGxxx'; None = all) and compute F/P/R.
    Summation is in RefOG order (identical to upstream for the full set)."""
    ids = None if ids is None else set(ids)
    tp = fp = fn = 0.0
    n = ex = 0
    for r in per_refog:
        if ids is None or r["refog"] in ids:
            tp += r["TP"]; fp += r["FP"]; fn += r["FN"]; n += 1; ex += int(r["exact"])
    F, P, R = prf(tp, fp, fn)
    return dict(F=F, precision=P, recall=R, F_frac=F / 100, precision_frac=P / 100, recall_frac=R / 100,
                TP=tp, FP=fp, FN=fn, n_refogs=n, exact_refogs=ex)


def score_orthogroups(bench_dir, ogs_file, basic=False, split_json=DEFAULT_SPLIT, verify=True, quiet=False):
    ob = load_official(bench_dir)
    d_ref = os.path.join(bench_dir, "RefOGs") + os.sep
    ref = ob.read_refogs(d_ref)
    unc = ob.read_uncertain_refogs(d_ref + "low_certainty_assignments/")
    exp = ob.get_expected_genes()
    pred = ob.read_orthogroups(ogs_file, exp, ob.get_n_col_skip(ogs_file)) if basic else ob.read_orthogroups_smart(ogs_file)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ob.check_orthogroups(pred, exp)
    if not quiet and buf.getvalue().strip():
        sys.stderr.write(buf.getvalue())
    per = per_refog_counts(ref, unc, pred)
    res = dict(per_refog=per)
    res["all"] = subset_scores(per)
    if verify:
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                f, p, r = ob.calculate_benchmarks_pairwise(ref, unc, pred)
            n_exact = int(re.search(r"(\d+) Orthogroups exactly correct", buf.getvalue()).group(1))
        except ZeroDivisionError:
            f = p = r = 0.0; n_exact = 0
        a = res["all"]
        ok = (f == a["F"] and p == a["precision"] and r == a["recall"]) and n_exact == a["exact_refogs"]
        res["official"] = dict(F=f, precision=p, recall=r, exact_refogs=n_exact)
        res["verify"] = dict(identical_to_official=bool(ok),
                             max_abs_diff=max(abs(f - a["F"]), abs(p - a["precision"]), abs(r - a["recall"])))
        if not ok and res["verify"]["max_abs_diff"] > 1e-9:
            raise AssertionError("per-RefOG sums differ from official benchmark.py: %r vs %r" % (res["official"], a))
    if split_json and os.path.exists(split_json):
        sp = json.load(open(split_json))
        res["split_file"] = os.path.abspath(split_json)
        res["dev"] = subset_scores(per, sp["dev"])
        res["test"] = subset_scores(per, sp["test"])
        assert res["dev"]["n_refogs"] + res["test"]["n_refogs"] == 70, "dev+test must cover the 70 RefOGs"
    res["n_pred_groups"] = len(pred)
    res["n_genes_assigned"] = sum(len(x) for x in pred)
    return res


def bootstrap_delta_f(perA, perB, ids=None, B=2000, seed=11):
    """Paired bootstrap over RefOGs of dF = F(A) - F(B), in percent points, on subset `ids`.
    perA/perB: per_refog lists (same RefOG order). RefOGs are resampled with replacement (same indices for both
    methods); F of each resample = F computed from the summed per-RefOG TP/FP/FN. numpy default_rng(seed).
    Returns dict(delta_point, mean, se, ci95=[lo,hi] (2.5/97.5 percentiles), p_le0 (share of resamples with
    dF<=0), B, seed, n_refogs, F_A, F_B)."""
    import numpy as np
    ids = None if ids is None else set(ids)
    idx = [i for i, r in enumerate(perA) if ids is None or r["refog"] in ids]
    assert [perA[i]["refog"] for i in idx] == [perB[i]["refog"] for i in idx]
    def arr(per):
        return np.array([[per[i]["TP"], per[i]["FP"], per[i]["FN"]] for i in idx])
    a, b = arr(perA), arr(perB)
    n = len(idx)
    rng = np.random.default_rng(seed)
    S = rng.integers(0, n, size=(B, n))
    def f_of(x):  # x: (B,3) sums
        tp, fp, fn = x[:, 0], x[:, 1], x[:, 2]
        with np.errstate(divide="ignore", invalid="ignore"):
            p = tp / (tp + fp); r = tp / (tp + fn)
            f = 100. * (2 * p * r / (p + r))
        return np.where(tp > 0, f, 0.0)
    fa = f_of(a[S].sum(axis=1)); fb = f_of(b[S].sum(axis=1))
    d = fa - fb
    fA = prf(*a.sum(axis=0))[0]; fB = prf(*b.sum(axis=0))[0]
    lo, hi = np.percentile(d, [2.5, 97.5])
    return dict(delta_point=fA - fB, mean=float(d.mean()), se=float(d.std(ddof=1)), ci95=[float(lo), float(hi)],
                p_le0=float((d <= 0).mean()), B=B, seed=seed, n_refogs=n, F_A=fA, F_B=fB)


def ni_verdict(ci_lo, margin):
    """Non-inferiority verdict for dF = F(A) - F(B) (percentage points): comparable iff lower 95% bound >= -margin."""
    if ci_lo >= -margin:
        return "comparable (non-inferior%s)" % ("; CI > 0" if ci_lo > 0 else "")
    return "not shown non-inferior"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bench_dir", nargs="?"); p.add_argument("ogs", nargs="?")
    p.add_argument("--basic", action="store_true")
    p.add_argument("--json"); p.add_argument("--label", default=None)
    p.add_argument("--split", default=DEFAULT_SPLIT)
    p.add_argument("--per-refog-tsv")
    p.add_argument("--no-verify", action="store_true")
    p.add_argument("--compare", nargs=2, metavar=("A.json", "B.json"))
    p.add_argument("--subset", default="all", choices=["all", "dev", "test"])
    p.add_argument("-B", type=int, default=2000); p.add_argument("--seed", type=int, default=11)
    p.add_argument("--ni-margin", type=float, default=None,
                   help="non-inferiority margin in percentage points (adds 'ni_verdict' to --compare output)")
    a = p.parse_args()
    if a.compare:
        A, Bj = (json.load(open(x)) for x in a.compare)
        ids = None
        if a.subset != "all":
            ids = json.load(open(a.split))[a.subset]
        r = bootstrap_delta_f(A["per_refog"], Bj["per_refog"], ids, a.B, a.seed)
        r.update(file_A=a.compare[0], file_B=a.compare[1], subset=a.subset)
        if a.ni_margin is not None:
            r.update(ni_margin_pp=a.ni_margin, ni_verdict=ni_verdict(r["ci95"][0], a.ni_margin))
        print(json.dumps(r, indent=1))
        return
    if not (a.bench_dir and a.ogs):
        p.error("bench_dir and ogs required")
    res = score_orthogroups(a.bench_dir, a.ogs, a.basic, a.split, verify=not a.no_verify)
    res["label"] = a.label or os.path.basename(a.ogs)
    al = res["all"]
    # backward compatible top-level keys (rounded as before)
    res.update(F=round(al["F"], 2), precision=round(al["precision"], 2), recall=round(al["recall"], 2),
               exact_refogs=al["exact_refogs"])
    print("label\tsubset\tF\tprecision\trecall\texact_refogs\tTP\tFP\tFN")
    for s in ("all", "dev", "test"):
        if s in res:
            x = res[s]
            print("%s\t%s\t%.2f\t%.2f\t%.2f\t%d/%d\t%.1f\t%.1f\t%.1f" % (res["label"], s, x["F"], x["precision"],
                  x["recall"], x["exact_refogs"], x["n_refogs"], x["TP"], x["FP"], x["FN"]))
    if "verify" in res:
        print("verify: per-RefOG sum identical to official benchmark.py: %s (max |diff| %.3g)" %
              (res["verify"]["identical_to_official"], res["verify"]["max_abs_diff"]))
    if a.json:
        json.dump(res, open(a.json, "w"), indent=1)
    if a.per_refog_tsv:
        with open(a.per_refog_tsv, "w") as fh:
            fh.write("refog\tn_genes\tTP\tFP\tFN\tn_pred_fragments\texact\n")
            for r in res["per_refog"]:
                fh.write("%s\t%d\t%.6f\t%.6f\t%.6f\t%d\t%d\n" % (r["refog"], r["n_genes"], r["TP"], r["FP"], r["FN"],
                                                              r["n_pred_fragments"], r["exact"]))


if __name__ == "__main__":
    main()
