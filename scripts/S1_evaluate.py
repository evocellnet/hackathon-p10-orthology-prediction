#!/usr/bin/env python
"""S1 evaluation: compare every guide-tree arm with (1) the NCBI taxdump tree and
(2) a published reference phylogeny -- EVALUATION ONLY, never a method input.

* QfO: the SwissTree species tree shipped with the QfO reference data
  (data/qfo/reference_data/2022.1/lineage_tree.phyloxml, 'SwissTree speciestree for
  QfO 2018 dataset', covers all 78 species of 2022_02). It is converted to Newick
  under results/S1_speciestree/eval_only/ (a directory the method never reads).
* OrthoBench: the 12-species NCBI tree is fully resolved and identical to the
  published vertebrate/bilaterian phylogeny (Ecdysozoa vs Chordata; Euarchontoglires
  incl. rodents; Clupeocephala), so the NCBI comparison doubles as the published one.

Metrics (ancestra.speciestree.compare_trees): nRF, ref_split_recall, conflict_frac,
quartet distance (resolved-only and all), overall and per domain
(Bacteria / Archaea / Eukaryota; QfO only).
"""
import glob
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from lxml import etree  # noqa: E402

from ancestra import speciestree as st  # noqa: E402
from ancestra.newick import Node, read_newick  # noqa: E402

RES = os.path.join(ROOT, "results", "S1_speciestree")
EVAL_REF = os.path.join(ROOT, "data", "qfo", "reference_data", "2022.1", "lineage_tree.phyloxml")


def phyloxml_tree(path: str, taxid_to_leaf: dict) -> Node:
    ns = {"p": "http://www.phyloxml.org"}
    doc = etree.parse(path)
    top = doc.find(".//p:phylogeny/p:clade", ns)

    def build(el, parent):
        tx = el.find("p:taxonomy", ns)
        tid = tx.findtext("p:id", namespaces=ns) if tx is not None else None
        nm = tx.findtext("p:scientific_name", namespaces=ns) if tx is not None else ""
        kids = el.findall("p:clade", ns)
        nd = Node(taxid_to_leaf.get(int(tid), nm) if (tid and not kids) else (nm or ""))
        if parent is not None:
            parent.add_child(nd)
        for k in kids:
            build(k, nd)
        return nd
    root = build(top, None)
    return st.collapse_unary(root)


LABEL = {}


def summarize(tree, ref, domains, tag):
    out = {"all": st.compare_trees(tree, ref)}
    out["conflicting_splits"] = [[LABEL.get(x, x) for x in sp] for sp in st.conflicting_splits(tree, ref)]
    if domains:
        for d in ("Bacteria", "Archaea", "Eukaryota"):
            leaves = [s for s, x in domains.items() if x == d]
            out[d] = st.compare_trees(tree, ref, leaves=leaves)
        # domain monophyly in query tree
        cl = set(st.clades(tree))
        out["domain_monophyly"] = {d: frozenset(s for s, x in domains.items() if x == d) in cl
                                   for d in ("Bacteria", "Archaea", "Eukaryota")}
    return out


def main():
    metrics = {}
    os.makedirs(os.path.join(RES, "eval_only"), exist_ok=True)
    for ds in ("ob", "qfo"):
        d = os.path.join(RES, ds)
        if not os.path.exists(os.path.join(d, "ncbi.nwk")):
            continue
        sp = json.load(open(os.path.join(d, "species.json")))
        if ds == "qfo":   # readable labels for conflicting-split lists (NCBI names)
            _, _, nm, _ = st.load_taxdump(os.path.join(ROOT, "data", "taxonomy"),
                                          need_names_for=set(sp["leaf_taxid"].values()))
            LABEL.update({leaf: nm.get(t, leaf) for leaf, t in sp["leaf_taxid"].items()})
        else:
            LABEL.update({leaf: leaf.split(".")[0] for leaf in sp["leaf_taxid"]})
        domains = sp["domain"] if ds == "qfo" else None
        ncbi = read_newick(os.path.join(d, "ncbi.nwk"))
        refs = {"ncbi": ncbi}
        if ds == "qfo":
            lt = {v: k for k, v in sp["leaf_taxid"].items()}
            pub = phyloxml_tree(EVAL_REF, lt)
            missing = set(lt.values()) - set(pub.leaf_names())
            assert not missing, missing
            st.write_tree(pub, os.path.join(RES, "eval_only", "qfo_swisstree_lineage_tree.nwk"), lengths=False)
            refs["published"] = pub
        metrics[ds] = {"refs": {k: {"n_leaves": len(v.leaves()),
                                    "n_internal": sum(1 for n in v.preorder() if n.children)}
                                for k, v in refs.items()}, "trees": {}}
        if ds == "qfo":
            metrics[ds]["refs"]["ncbi_vs_published"] = summarize(ncbi, refs["published"], domains, "")
        for f in sorted(glob.glob(os.path.join(d, "*.nwk"))):
            name = os.path.basename(f)[:-4]
            if name == "ncbi":
                continue
            t = read_newick(f)
            metrics[ds]["trees"][name] = {rn: summarize(t, r, domains, name) for rn, r in refs.items()}
            print(ds, name, {rn: round(v["all"]["nrf"], 3) for rn, v in metrics[ds]["trees"][name].items()},
                  flush=True)
    with open(os.path.join(RES, "metrics.json"), "w") as fh:
        json.dump(metrics, fh, indent=1)




def table(metrics_path=os.path.join(RES, "metrics.json"), out_md=os.path.join(RES, "summary.md")):
    """Markdown/TSV summary generated from metrics.json (perturbed/random: mean over seeds)."""
    import re
    from collections import defaultdict
    import numpy as np
    m = json.load(open(metrics_path))
    lines = []
    rows_tsv = ["ds\ttree\tref\tscope\tnrf\tref_split_recall\tconflict_frac\tqd_resolved\tqd_all\tn_seeds"]
    for ds, dm in m.items():
        groups = defaultdict(list)
        for name in dm["trees"]:
            g = re.sub(r"_s(11|23|47)$", "_s*", name)
            groups[g].append(name)
        refs = [r for r in ("ncbi", "published") if r in dm["refs"]]
        scopes = ["all", "Bacteria", "Archaea", "Eukaryota"] if ds == "qfo" else ["all"]
        lines.append(f"\n### {ds}\n")
        hdr = "| tree | " + " | ".join(f"{r} nRF / recall / conflict / qd_res" for r in refs)
        if ds == "qfo":
            hdr += " | published nRF (conflict) Bact / Arch / Euk | domains monophyletic (B,A,E)"
        lines.append(hdr + " |")
        lines.append("|" + "---|" * (1 + len(refs) + (2 if ds == "qfo" else 0)))
        for g, names in groups.items():
            cells = [g]
            for r in refs:
                vals = {k: np.mean([dm["trees"][n][r]["all"][k] for n in names]) for k in
                        ("nrf", "ref_split_recall", "conflict_frac")}
                qd = np.mean([dm["trees"][n][r]["all"]["quartet"]["qd_resolved"] for n in names])
                cells.append(f"{vals['nrf']:.3f} / {vals['ref_split_recall']:.3f} / "
                             f"{vals['conflict_frac']:.3f} / {qd:.3f}")
            for r in refs:
                for sc in scopes:
                    for n in names:
                        x = dm["trees"][n][r].get(sc, {})
                        if "nrf" in x:
                            rows_tsv.append(f"{ds}\t{n}\t{r}\t{sc}\t{x['nrf']:.4f}\t{x['ref_split_recall']:.4f}\t"
                                            f"{x['conflict_frac']:.4f}\t{x['quartet']['qd_resolved']:.4f}\t"
                                            f"{x['quartet']['qd_all']:.4f}\t{len(names)}")
            if ds == "qfo":
                r = "published"
                cells.append(" / ".join(f"{np.mean([dm['trees'][n][r][sc]['nrf'] for n in names]):.2f}"
                                        f" ({np.mean([dm['trees'][n][r][sc]['conflict_frac'] for n in names]):.2f})"
                                        for sc in ("Bacteria", "Archaea", "Eukaryota")))
                mono = [np.mean([dm["trees"][n]["ncbi"]["domain_monophyly"][d] for n in names])
                        for d in ("Bacteria", "Archaea", "Eukaryota")]
                cells.append(",".join(f"{x:.2g}" for x in mono))
            lines.append("| " + " | ".join(cells) + " |")
        if ds == "qfo" and "ncbi_vs_published" in dm["refs"]:
            x = dm["refs"]["ncbi_vs_published"]
            lines.append(f"\nNCBI vs published (SwissTree lineage tree): nRF {x['all']['nrf']:.3f}, "
                         f"recall {x['all']['ref_split_recall']:.3f}, conflict {x['all']['conflict_frac']:.3f}; "
                         f"splits NCBI {x['all']['n_splits_query']}, published {x['all']['n_splits_ref']}, "
                         f"binary max {x['all']['max_splits']}.")
    open(out_md, "w").write("\n".join(lines) + "\n")
    open(out_md.replace(".md", ".tsv"), "w").write("\n".join(rows_tsv) + "\n")
    print("\n".join(lines))


def manifest():
    """results/S1_speciestree/manifest.json: inputs (checksums), commands, seeds, runtimes, outputs."""
    import hashlib
    import subprocess
    import time

    def sha(p):
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for b in iter(lambda: fh.read(1 << 20), b""):
                h.update(b)
        return h.hexdigest()
    tax = os.path.join(ROOT, "data", "taxonomy")
    try:
        commit = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True,
                                text=True).stdout.strip()
    except Exception:  # noqa: BLE001
        commit = None
    man = {
        "subtask": "S1 guide species trees", "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "git_commit_at_creation": commit, "seeds": [11, 23, 47],
        "code": {p: sha(os.path.join(ROOT, p)) for p in
                 ("ancestra/speciestree.py", "scripts/S1_run.py", "scripts/S1_evaluate.py",
                  "tests/test_speciestree.py")},
        "inputs": {
            "ncbi_taxdump": {"url": "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
                             "md5": open(os.path.join(tax, "taxdump.tar.gz.md5")).read().split()[0],
                             "md5_verified_local": None,
                             "downloaded": open(os.path.join(tax, "download_date.txt")).read().strip()},
            "qfo_proteomes": "data/qfo/proteomes_2022_02_updZEBRAFISH (md5sums.txt there; DATA.md 3)",
            "orthobench_proteomes": "data/orthobench/Open_Orthobench/BENCHMARKS/Input (DATA.md 1)",
            "s0_families_ob": {"path": "work/S0/2sA_dc30_cc_sh30_RA/families.tsv",
                               "sha256": sha(os.path.join(ROOT, "work/S0/2sA_dc30_cc_sh30_RA/families.tsv"))},
            "evaluation_only": {"qfo_published_tree": "data/qfo/reference_data/2022.1/lineage_tree.phyloxml",
                                "sha256": sha(EVAL_REF)},
        },
        "commands": [
            "python scripts/S1_run.py {ob,qfo} ncbi search armA armB perturb --threads 12",
            "python scripts/S1_run.py ob armB --threads 2 --max-copy 1 --max-mean-copy 1.0 --tag _sc1",
            "python scripts/S1_evaluate.py",
        ],
        "runtime": {}, "outputs": {},
    }
    md5 = hashlib.md5()
    with open(os.path.join(tax, "taxdump.tar.gz"), "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            md5.update(b)
    man["inputs"]["ncbi_taxdump"]["md5_verified_local"] = md5.hexdigest()
    for ds in ("ob", "qfo"):
        rp = os.path.join(RES, ds, "runtime.json")
        if os.path.exists(rp):
            man["runtime"][ds] = json.load(open(rp))
        for f in sorted(glob.glob(os.path.join(RES, ds, "*"))):
            man["outputs"][os.path.relpath(f, RES)] = sha(f)
    for f in ("metrics.json", "summary.md", "summary.tsv", "eval_only/qfo_swisstree_lineage_tree.nwk"):
        if os.path.exists(os.path.join(RES, f)):
            man["outputs"][f] = sha(os.path.join(RES, f))
    with open(os.path.join(RES, "manifest.json"), "w") as fh:
        json.dump(man, fh, indent=1)


if __name__ == "__main__":
    if len(sys.argv) == 1 or sys.argv[1] != "table":
        main()
    table()
    manifest()
