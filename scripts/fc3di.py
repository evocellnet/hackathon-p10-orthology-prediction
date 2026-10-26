#!/usr/bin/env python
"""Helpers for building 3Di strings from Foldcomp databases via `foldseek createdb --input-format 5`.

Subcommands
  select   <lookup> <acc_json|acc_txt> <keys_out> <map_out> [--ver v6]
           Stream a foldcomp .lookup (key \t name \t x), find AF-<acc>-F1-model_<ver>[.pdb] for each wanted accession,
           write the numeric keys (for `foldseek createsubdb --id-mode 0`) and a TSV key->acc.
  extract  <fs_db_prefix> <out_tsv>
           Read a foldseek DB produced by createdb (AA db + _ss) and write TSV: acc, aa, 3di  (F1 entries only)
  qfo      <fs_db_prefix> <qfo_fasta_dir> <out_fa> <out_stats_tsv> [--ref-3di <fasta>]
           QfO: keep entries whose AA == QfO sequence (exact); output >UPID|acc 3Di
  ob       <fs_db_prefix> <orthobench_input_dir> <uniprot_seq_tsv> <out_fa> <out_stats_tsv>
           OrthoBench: Ensembl protein -> exact AA match to an AFDB entry (via fs DB AA, plus UniProt seq->acc table)
"""
import sys, os, re, glob, json, argparse, mmap
import numpy as np


def read_mmdb(prefix):
    """Read an MMseqs/Foldseek DB of text entries -> dict key(int) -> str (trailing \\n\\0 stripped)."""
    idx = np.loadtxt(prefix + '.index', dtype=np.int64, ndmin=2) if os.path.getsize(prefix + '.index') else np.zeros((0, 3), np.int64)
    out = {}
    with open(prefix, 'rb') as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        for k, off, ln in idx:
            out[int(k)] = mm[off:off + ln - 2].decode()  # drop \n\0
    return out


def read_fs_db(prefix):
    """-> list of (acc, fragment, aa, ss) from a foldseek createdb output built from Foldcomp/AFDB entries."""
    aa = read_mmdb(prefix)
    ss = read_mmdb(prefix + '_ss')
    hd = read_mmdb(prefix + '_h')
    res = []
    for k, s in aa.items():
        m = re.match(r'AF-(.+)-F(\d+)-model_v\d+', hd[k])
        if m:
            res.append((m.group(1), int(m.group(2)), s, ss[k]))
    return res


def iter_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for l in f:
            if l.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name, seq = l[1:].strip().split()[0], []
            else:
                seq.append(l.strip())
    if name is not None:
        yield name, ''.join(seq)


def cmd_select(a):
    wanted = set()
    p = a.acc
    if p.endswith('.json'):
        d = json.load(open(p))
        for v in d.values():
            wanted.update(v)
    else:
        wanted.update(x.strip() for x in open(p) if x.strip())
    n = 0
    with open(a.lookup) as f, open(a.keys_out, 'w') as ko, open(a.map_out, 'w') as mo:
        for l in f:
            t = l.split('\t', 2)
            p = t[1].split('-')  # AF-<acc>-F<n>-model_vX[.pdb]
            if len(p) == 4 and p[2] == 'F1' and p[1] in wanted and (a.ver is None or p[3].startswith('model_' + a.ver)):
                ko.write(t[0] + '\n')
                mo.write(f'{t[0]}\t{p[1]}\n')
                n += 1
    print(f'wanted accessions {len(wanted)}; found F1 entries {n}', file=sys.stderr)


def cmd_extract(a):
    with open(a.out, 'w') as o:
        for acc, fr, s, ss in read_fs_db(a.prefix):
            if fr == 1:
                o.write(f'{acc}\t{s}\t{ss}\n')


def cmd_qfo(a):
    fs = {acc: (s, ss) for acc, fr, s, ss in read_fs_db(a.prefix) if fr == 1}
    ref = {}
    if a.ref_3di:
        for n, s in iter_fasta(a.ref_3di):
            ref[n] = s
    stats = []
    cmp_n = cmp_same = 0
    alias = {}
    if a.alias:
        for l in open(a.alias):
            x, p = l.split()
            alias.setdefault(x, []).append(p)
    src = open(a.out_fa + '.src.tsv', 'w')
    with open(a.out_fa, 'w') as o:
        for f in sorted(glob.glob(os.path.join(a.qfo_dir, '*.fasta'))):
            pid = os.path.basename(f).split('_')[0]
            tot = inaf = exact = 0
            for r, seq in iter_fasta(f):
                acc = r.split('|')[1]
                tot += 1
                hit = None
                for c in [acc] + alias.get(acc, []):
                    if c in fs:
                        inaf += 1 if c == acc else 0
                        if fs[c][0] == seq:
                            hit = c
                            break
                if hit:
                    if True:
                        exact += 1
                        o.write(f'>{pid}|{acc}\n{fs[hit][1]}\n')
                        src.write(f'{pid}|{acc}\t{hit}\t{"own_acc" if hit == acc else "merged_acc"}\n')
                        if f'{pid}|{acc}' in ref:
                            cmp_n += 1
                            cmp_same += ref[f'{pid}|{acc}'] == fs[hit][1]
            stats.append((pid, tot, inaf, exact))
    src.close()
    with open(a.out_stats, 'w') as o:
        o.write('proteome\tn_proteins\tn_in_foldcomp\tn_exact_seq\n')
        for s in stats:
            o.write('\t'.join(map(str, s)) + '\n')
    print('QfO: total', sum(s[1] for s in stats), 'in_foldcomp', sum(s[2] for s in stats), 'exact', sum(s[3] for s in stats))
    if a.ref_3di:
        print(f'3Di vs reference (Foldseek AFDB Proteome DB, overlap with exact AA): {cmp_same}/{cmp_n} identical')


def seq_norm(s):
    return s.strip().rstrip('*').upper()


def cmd_obcand(a):
    """OrthoBench Ensembl protein -> candidate UniProt accessions by exact sequence (UniProt taxon FASTAs + QfO FASTAs)."""
    import gzip
    ob = {}   # seq -> [(species, ensembl id)]
    for f in sorted(glob.glob(os.path.join(a.ob_dir, '*.fa'))):
        sp = os.path.basename(f).split('.')[0]
        for n, s in iter_fasta(f):
            ob.setdefault(seq_norm(s), []).append((sp, n))
    cand = {}  # seq -> set(acc)
    def scan(path):
        op = gzip.open if path.endswith('.gz') else open
        name, seq = None, []
        def flush():
            if name is not None:
                s = ''.join(seq)
                if s in ob:
                    cand.setdefault(s, set()).add(name.split('|')[1] if '|' in name else name)
        with op(path, 'rt') as f:
            for l in f:
                if l.startswith('>'):
                    flush(); name, seq = l[1:].strip().split()[0], []
                else:
                    seq.append(l.strip())
            flush()
    for p in a.fasta:
        scan(p)
    accs = set()
    with open(a.out_map, 'w') as o:
        for s, ids in ob.items():
            c = sorted(cand.get(s, []))
            accs.update(c)
            for sp, n in ids:
                o.write(f'{sp}\t{n}\t{",".join(c)}\n')
    with open(a.out_acc, 'w') as o:
        o.write(''.join(x + '\n' for x in sorted(accs)))
    print('OrthoBench proteins', sum(len(v) for v in ob.values()), 'with >=1 UniProt candidate',
          sum(len(v) for s, v in ob.items() if s in cand), 'candidate accessions', len(accs), file=sys.stderr)


def cmd_ob(a):
    fs = {acc: (s, ss) for acc, fr, s, ss in read_fs_db(a.prefix) if fr == 1}
    ref = {}
    if a.ref_3di:
        for n, s in iter_fasta(a.ref_3di):
            ref[n.split('|')[0]] = s
    stats = {}
    cmp_n = cmp_same = 0
    seqs = {}
    for f in sorted(glob.glob(os.path.join(a.ob_dir, '*.fa'))):
        for n, s in iter_fasta(f):
            seqs[n] = seq_norm(s)
    src = open(a.out_fa + '.src.tsv', 'w')
    with open(a.out_fa, 'w') as o:
        for l in open(a.map):
            sp, n, c = l.rstrip('\n').split('\t')
            st = stats.setdefault(sp, [0, 0, 0, 0])
            st[0] += 1
            accs = [x for x in c.split(',') if x]
            if accs:
                st[1] += 1
            hit = None
            for x in accs:
                if x in fs and fs[x][0] == seqs[n]:
                    hit = x; break
            if hit is None and any(x in fs for x in accs):
                st[2] += 1   # in DB but sequence differs
            if hit:
                st[3] += 1
                o.write(f'>{n}|{hit}\n{fs[hit][1]}\n')
                src.write(f'{n}|{hit}\t{hit}\tseq\n')
                if n in ref:
                    cmp_n += 1; cmp_same += ref[n] == fs[hit][1]
    src.close()
    with open(a.out_stats, 'w') as o:
        o.write('species\tn_proteins\tn_uniprot_seq_match\tn_in_db_but_seq_differs\tn_exact_3di\n')
        for sp, st in sorted(stats.items()):
            o.write(sp + '\t' + '\t'.join(map(str, st)) + '\n')
    print('OrthoBench: total', sum(s[0] for s in stats.values()), 'exact', sum(s[3] for s in stats.values()))
    if a.ref_3di:
        print(f'3Di vs previous exact-match file: {cmp_same}/{cmp_n} identical')


def _scan_range(args):
    path, start, end, wanted = args
    hits = []
    with open(path, 'rb') as f:
        # align start to the beginning of an entry (after a \0)
        if start > 0:
            f.seek(start - 1)
            buf = b''
            pos = start - 1
            while True:
                chunk = f.read(1 << 16)
                i = chunk.find(b'\0')
                if i >= 0:
                    start = pos + i + 1
                    break
                pos += len(chunk)
                if not chunk:
                    return hits
        f.seek(start)
        pos = start
        rest = b''
        restpos = start
        while pos < end or rest:
            chunk = f.read(64 << 20)
            if not chunk and not rest:
                break
            data = rest + chunk
            base = restpos
            last = data.rfind(b'\0')
            if last < 0:
                rest, restpos = data, base
                if not chunk:
                    break
                continue
            body, rest = data[:last + 1], data[last + 1:]
            restpos = base + last + 1
            off = base
            for e in body.split(b'\0'):
                if not e:
                    continue
                s_ = e[:-1] if e.endswith(b'\n') else e
                if s_ in wanted:
                    hits.append((off, s_.decode()))
                off += len(e) + 1
                if off >= end:
                    break
            pos = off
            if off >= end:
                break
            if not chunk:
                break
    return hits


def cmd_seqscan(a):
    """Find AFDB entries whose AA sequence equals one of the wanted sequences. a.db = foldseek AA db data file."""
    import multiprocessing as mp
    wanted = set()
    for p in a.fasta:
        for n, s in iter_fasta(p):
            wanted.add(seq_norm(s).encode())
    size = os.path.getsize(a.db)
    n = a.chunks
    bounds = [size * i // n for i in range(n + 1)]
    with mp.Pool(a.procs) as pool:
        res = pool.map(_scan_range, [(a.db, bounds[i], bounds[i + 1], wanted) for i in range(n)], chunksize=1)
    hits = {}
    for r in res:
        for off, s in r:
            hits[off] = s
    # offsets -> keys via index
    offs = set(hits)
    out = []
    with open(a.db + '.index') as f:
        for l in f:
            k, o, ln = l.split()
            if int(o) in offs:
                out.append((k, hits[int(o)]))
    with open(a.out, 'w') as o:
        for k, s in out:
            o.write(f'{k}\t{s}\n')
    print(f'wanted {len(wanted)} sequences; AFDB entries with identical sequence: {len(out)}', file=sys.stderr)


def cmd_wanted(a):
    """Union of wanted accessions: all QfO accessions + OrthoBench candidate accessions + primary accessions of QfO
    accessions that are now secondary (UniProt sec_ac.txt). Writes wanted_acc.txt and qfo_alias.tsv (qfo_acc -> primary)."""
    qfo = json.load(open(a.qfo_json))
    allq = {x for v in qfo.values() for x in v}
    sec = {}
    started = False
    for l in open(a.sec_ac):
        if l.startswith('Secondary AC'):
            started = True
            continue
        p = l.split()
        if started and len(p) == 2 and re.fullmatch(r'[A-Z][A-Z0-9]{5,9}', p[0]):
            sec.setdefault(p[0], []).append(p[1])
    wanted = set(allq)
    n = 0
    with open(a.out_alias, 'w') as o:
        for x in sorted(allq):
            for pr in sec.get(x, []):
                o.write(f'{x}\t{pr}\n'); wanted.add(pr); n += 1
    for p in a.extra:
        wanted.update(l.strip() for l in open(p) if l.strip())
    with open(a.out, 'w') as o:
        o.write(''.join(x + '\n' for x in sorted(wanted)))
    print('wanted accessions', len(wanted), 'qfo aliases', n, file=sys.stderr)


def read_by_keys(prefix, keys):
    """Read selected keys from an MMseqs-style DB (index scan, then direct reads)."""
    want = set(keys)
    ent = {}
    with open(prefix + '.index') as f:
        for l in f:
            k, o, ln = l.split()
            if k in want:
                ent[k] = (int(o), int(ln))
    out = {}
    with open(prefix, 'rb') as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        for k, (o, ln) in ent.items():
            out[k] = mm[o:o + ln - 2].decode()
    return out


def cmd_fsexact(a):
    """Exact-sequence route on a Foldseek AFDB DB (afdb/afdb_ss/lookup): find AFDB entries whose sequence equals a QfO or
    OrthoBench protein (any accession), take their true 3Di. Prefers the protein's own UniProt accession."""
    import subprocess
    qfo_files = sorted(glob.glob(os.path.join(a.qfo_dir, '*.fasta')))
    ob_files = sorted(glob.glob(os.path.join(a.ob_dir, '*.fa')))
    hits_tsv = a.out_prefix + '.seqhits.tsv'
    if not os.path.exists(hits_tsv):
        ns = argparse.Namespace(db=a.db, out=hits_tsv, fasta=qfo_files + ob_files, procs=a.procs, chunks=a.chunks)
        cmd_seqscan(ns)
    import pickle
    pk = a.out_prefix + '.byseq.pkl'
    if os.path.exists(pk):
        byseq = pickle.load(open(pk, 'rb'))
    else:
        keys = {}
        for l in open(hits_tsv):
            k, s = l.rstrip('\n').split('\t'); keys[k] = s
        names = {}
        with open(a.db + '.lookup') as f:
            for l in f:
                t = l.split('\t', 2)
                if t[0] in keys:
                    names[t[0]] = t[1]
        ss = read_by_keys(a.db + '_ss', keys)
        byseq = {}
        for k, s in keys.items():
            m = re.match(r'AF-(.+)-F(\d+)-model', names[k])
            if m:
                byseq.setdefault(s, []).append((int(m.group(2)), m.group(1), ss[k]))
        for s in byseq:
            byseq[s].sort()
        pickle.dump(byseq, open(pk, 'wb'))
    obcand = {}
    if a.ob_cand:
        for l in open(a.ob_cand):
            sp, n, c = l.rstrip('\n').split('\t')
            obcand[n] = set(x for x in c.split(',') if x)
    stat_q, stat_o = {}, {}
    with open(a.out_prefix + '.qfo.fa', 'w') as o, open(a.out_prefix + '.qfo.src.tsv', 'w') as so:
        for f in qfo_files:
            pid = os.path.basename(f).split('_')[0]
            st = stat_q.setdefault(pid, [0, 0, 0])
            for r, seq in iter_fasta(f):
                acc = r.split('|')[1]
                st[0] += 1
                h = byseq.get(seq_norm(seq))
                if not h:
                    continue
                own = [x for x in h if x[1] == acc]
                fr, hacc, d = (own or h)[0]
                st[1 if own else 2] += 1
                o.write(f'>{pid}|{acc}\n{d}\n')
                so.write(f'{pid}|{acc}\t{hacc}\t{"own_acc" if own else "seq_other_acc"}\n')
    with open(a.out_prefix + '.ob.fa', 'w') as o, open(a.out_prefix + '.ob.src.tsv', 'w') as so:
        for f in ob_files:
            sp = os.path.basename(f).split('.')[0]
            st = stat_o.setdefault(sp, [0, 0])
            for n, seq in iter_fasta(f):
                st[0] += 1
                h = byseq.get(seq_norm(seq))
                if not h:
                    continue
                pref = [x for x in h if x[1] in obcand.get(n, ())]
                fr, hacc, d = (pref or h)[0]
                st[1] += 1
                o.write(f'>{n}|{hacc}\n{d}\n')
                so.write(f'{n}|{hacc}\t{hacc}\t{"uniprot_candidate" if pref else "seq_other_acc"}\n')
    with open(a.out_prefix + '.qfo.stats.tsv', 'w') as o:
        o.write('proteome\tn_proteins\tn_exact_own_acc\tn_exact_other_acc\n')
        for k, v in sorted(stat_q.items()):
            o.write(k + '\t' + '\t'.join(map(str, v)) + '\n')
    with open(a.out_prefix + '.ob.stats.tsv', 'w') as o:
        o.write('species\tn_proteins\tn_exact_seq\n')
        for k, v in sorted(stat_o.items()):
            o.write(k + '\t' + '\t'.join(map(str, v)) + '\n')
    print('fsexact QfO', sum(v[0] for v in stat_q.values()), 'exact', sum(v[1] + v[2] for v in stat_q.values()),
          '| OB', sum(v[0] for v in stat_o.values()), 'exact', sum(v[1] for v in stat_o.values()))


def cmd_merge(a):
    """Merge several 3Di FASTAs by priority (first label wins). Writes <out>.fa, <out>.source.tsv, <out>.coverage.tsv."""
    srcs = []
    for x in a.src:
        lab, p = x.split('=', 1)
        if os.path.exists(p):
            srcs.append((lab, p))
        else:
            print('WARNING: missing source', lab, p, file=sys.stderr)
    # universe: id -> group
    uni = {}
    if a.kind == 'qfo':
        for f in sorted(glob.glob(os.path.join(a.dir, '*.fasta'))):
            g = os.path.basename(f).rsplit('.', 1)[0]
            for r, _ in iter_fasta(f):
                uni[g.split('_')[0] + '|' + r.split('|')[1]] = g
    else:
        for f in sorted(glob.glob(os.path.join(a.dir, '*.fa'))):
            g = os.path.basename(f).split('.')[0]
            for n, _ in iter_fasta(f):
                uni[n] = g
    keyof = (lambda i: i) if a.kind == 'qfo' else (lambda i: i.split('|')[0])
    done = {}   # key -> (label, id, 3di, srcinfo)
    for lab, p in srcs:
        info = {}
        sp = p + '.src.tsv' if os.path.exists(p + '.src.tsv') else p.replace('.fa', '.src.tsv')
        if os.path.exists(sp):
            for l in open(sp):
                t = l.rstrip('\n').split('\t')
                info[t[0]] = t[1:]
        for i, s in iter_fasta(p):
            k = keyof(i)
            if k in uni and k not in done:
                inf = info.get(i, ['', ''])
                done[k] = (lab, i, s, inf)
    cnt = {}
    with open(a.out + '.fa', 'w') as o, open(a.out + '.source.tsv', 'w') as so:
        so.write('id\tsource\tmatched_accession\tmatch_kind\n')
        for k, (lab, i, s, inf) in done.items():
            o.write(f'>{i}\n{s}\n')
            so.write(f'{i}\t{lab}\t{inf[0] if inf else ""}\t{inf[1] if len(inf) > 1 else ""}\n')
            c = cnt.setdefault(uni[k], {})
            c[lab] = c.get(lab, 0) + 1
    tot = {}
    for k, g in uni.items():
        tot[g] = tot.get(g, 0) + 1
    labs = [l for l, _ in srcs]
    with open(a.out + '.coverage.tsv', 'w') as o:
        o.write('group\tn_proteins\tn_with_3di\tfrac\t' + '\t'.join('n_' + l for l in labs) + '\n')
        for g in sorted(tot):
            c = cnt.get(g, {})
            n = sum(c.values())
            o.write(f'{g}\t{tot[g]}\t{n}\t{n/tot[g]:.4f}\t' + '\t'.join(str(c.get(l, 0)) for l in labs) + '\n')
        T = sum(tot.values()); N = len(done)
        o.write(f'TOTAL\t{T}\t{N}\t{N/T:.4f}\t' + '\t'.join(str(sum(c.get(l, 0) for c in cnt.values())) for l in labs) + '\n')
    print(a.kind, 'total', sum(tot.values()), 'with 3Di', len(done), {l: sum(c.get(l, 0) for c in cnt.values()) for l in labs})


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    s = sp.add_parser('select'); s.add_argument('lookup'); s.add_argument('acc'); s.add_argument('keys_out'); s.add_argument('map_out'); s.add_argument('--ver')
    s = sp.add_parser('extract'); s.add_argument('prefix'); s.add_argument('out')
    s = sp.add_parser('qfo'); s.add_argument('prefix'); s.add_argument('qfo_dir'); s.add_argument('out_fa'); s.add_argument('out_stats'); s.add_argument('--ref-3di'); s.add_argument('--alias')
    s = sp.add_parser('obcand'); s.add_argument('ob_dir'); s.add_argument('out_map'); s.add_argument('out_acc'); s.add_argument('fasta', nargs='+')
    s = sp.add_parser('ob'); s.add_argument('prefix'); s.add_argument('ob_dir'); s.add_argument('map'); s.add_argument('out_fa'); s.add_argument('out_stats'); s.add_argument('--ref-3di')
    s = sp.add_parser('seqscan'); s.add_argument('db'); s.add_argument('out'); s.add_argument('fasta', nargs='+'); s.add_argument('--procs', type=int, default=16); s.add_argument('--chunks', type=int, default=64)
    s = sp.add_parser('wanted'); s.add_argument('qfo_json'); s.add_argument('sec_ac'); s.add_argument('out'); s.add_argument('out_alias'); s.add_argument('extra', nargs='*')
    s = sp.add_parser('fsexact'); s.add_argument('db'); s.add_argument('qfo_dir'); s.add_argument('ob_dir'); s.add_argument('out_prefix'); s.add_argument('--procs', type=int, default=16); s.add_argument('--chunks', type=int, default=64)
    s = sp.add_parser('merge'); s.add_argument('kind', choices=['qfo', 'ob']); s.add_argument('dir'); s.add_argument('out'); s.add_argument('src', nargs='+')
    a = ap.parse_args()
    {'select': cmd_select, 'extract': cmd_extract, 'qfo': cmd_qfo, 'obcand': cmd_obcand, 'ob': cmd_ob, 'seqscan': cmd_seqscan, 'wanted': cmd_wanted, 'fsexact': cmd_fsexact, 'merge': cmd_merge}[a.cmd](a)
