#!/usr/bin/env python
"""L2-gpu task 2 micro-benchmark on real local-tree jobs dumped from the frozen OrthoBench accurate merge
(first two FastTree calls of the constrained pass): equality and speed of
  (1) FAMSA CLI vs in-process pyfamsa alignments,
  (2) one FastTree process per job vs FastTree -n over batches (same alignments),
single-threaded, sequential (CPU time = wall).
    python scripts/l2gpu_t2_bench.py PAYLOADS.pkl OUT.json [--max-jobs N]
"""
import json
import os
import pickle
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ancestra.localtree import _famsa_text, _fasta_to_phylip  # noqa: E402


def ft_single(aln_txt, td):
    p = os.path.join(td, "a.fa")
    open(p, "w").write(aln_txt)
    return subprocess.run(["FastTree", "-lg", "-quiet", "-nosupport", p], check=True, capture_output=True,
                          text=True).stdout.strip()


def ft_batch(alns, td):
    p = os.path.join(td, "all.phy")
    with open(p, "w") as fh:
        for a in alns:
            fh.write(_fasta_to_phylip(a))
    r = subprocess.run(["FastTree", "-lg", "-quiet", "-nosupport", "-n", str(len(alns)), p], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(r.stderr[-1500:])
    return [x.strip() + ";" for x in r.stdout.split(";") if x.strip()]


def topo(nwk):
    """Unrooted split set of a Newick string (leaf names only) for topology comparison."""
    from ancestra.localtree import newick_to_edges
    names = sorted(set(t for t in nwk.replace("(", ",").replace(")", ",").replace(";", "").split(",") if t)
                   if False else [])
    import re
    leaves = re.findall(r"(L\d+)", nwk)
    idx = {l: i for i, l in enumerate(sorted(set(leaves), key=lambda s: int(s[1:])))}
    edges, n_nodes = newick_to_edges(nwk, idx)
    adj = {}
    for a, b in edges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    m = len(idx)
    splits = set()
    for a, b in edges:
        seen, st, side = {a}, [b], set()
        seen.add(b)
        while st:
            x = st.pop()
            if x < m:
                side.add(x)
            for y in adj.get(x, []):
                if y not in seen:
                    seen.add(y)
                    st.append(y)
        s = frozenset(side)
        if 1 < len(s) < m - 1:
            splits.add(min(s, frozenset(range(m)) - s, key=lambda z: sorted(z)))
    return splits


def main():
    pp, out = sys.argv[1], sys.argv[2]
    mx = int(sys.argv[sys.argv.index("--max-jobs") + 1]) if "--max-jobs" in sys.argv else 10 ** 9
    calls = pickle.load(open(pp, "rb"))
    jobs = [j for c in calls for j in c][:mx]
    res = {"n_jobs": len(jobs)}
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(out))) as td:
        t0 = time.time()
        cli = [_famsa_text(n, s, td, False) for n, s in jobs]
        res["famsa_cli_s"] = round(time.time() - t0, 2)
        t0 = time.time()
        inp, bad = [], 0
        for n, s in jobs:
            try:
                inp.append(_famsa_text(n, s, td, True))
            except ValueError:
                inp.append(None)
                bad += 1
        res["famsa_inproc_s"] = round(time.time() - t0, 2)
        res["inproc_rejected"] = bad
        norm = lambda t: "".join(l.strip() if not l.startswith(">") else "\n" + l for l in t.splitlines())  # noqa
        res["inproc_identical_alignments"] = sum(1 for a, b in zip(cli, inp) if b is not None and norm(a) == norm(b))
        t0 = time.time()
        single = [ft_single(a, td) for a in cli]
        res["fasttree_single_s"] = round(time.time() - t0, 2)
        for bs in (50, 200):
            t0 = time.time()
            batched = []
            for i in range(0, len(cli), bs):
                batched += ft_batch(cli[i:i + bs], td)
            res[f"fasttree_batch{bs}_s"] = round(time.time() - t0, 2)
            res[f"batch{bs}_identical_newick"] = sum(a == b for a, b in zip(single, batched))
            res[f"batch{bs}_identical_topology"] = sum(topo(a) == topo(b) for a, b in zip(single, batched))
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
