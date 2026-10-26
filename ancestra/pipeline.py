"""End-to-end core pipeline (v0, sequence only): proteomes + guide tree -> HOGs.

    python -m ancestra run --proteomes DIR --tree tree.nwk [--families fam.tsv] \\
        --out OUTDIR [--threads 24] [--arm max] [--config cfg.json]

Stages (each checkpointed in ``OUTDIR/stages/<stage>/`` with a ``DONE`` file
holding a hash of the stage's config and of its upstream stages; outputs are
written to a temp dir and renamed, so a killed run never leaves a half stage;
a rerun skips every stage whose key is unchanged):

1. ``input``    gene table (integer IDs), species table, pruned guide tree
2. ``families`` gene -> family (given TSV, or S0 default clustering via
                ``ancestra.families`` if importable, else DIAMOND all-vs-all +
                connected components)
3. ``scores``   S2 within-family scores (``ancestra.scoring``)
4. ``merge``    S3 bottom-up merge (``ancestra.merge``)
5. ``output``   OUTDIR/hogs.orthoxml (OrthoXML 0.5, XSD-validated),
                OUTDIR/pairs.tsv (QfO flattening of the HOGs; secondary export),
                OUTDIR/orthobench_ogs.txt (root-level groups, OrthoBench format)

``OUTDIR/manifest.json``: git commit, full config, md5 of inputs, per-stage
wall time and peak RSS (own process VmHWM reset per stage; children = max RSS
of any finished child process so far), per-node merge statistics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pickle
import platform
import resource
import shutil
import subprocess
import sys
import time
from dataclasses import asdict

import numpy as np

from .io import (GeneTable, families_from_mapping, load_families, load_guide_tree,
                 load_proteomes, write_fasta, write_gene_table)
from .merge import ARMS, MergeConfig, MergeResult, run_merge, to_orthoxml
from .orthoxml import hogs_to_pairs, to_orthobench_ogs, validate_xsd, write_orthoxml
from .scoring import ScoreTable, ScoringConfig, score_with_rescue

__all__ = ["run_pipeline", "main", "DEFAULT_FAMILY_CONFIG"]

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# S0 default = `2sA_dc30_cc_sh30_RA` (docs/S0_families.md, selected on OrthoBench DEV RefOGs):
# DIAMOND cluster 30% id / 50% member cover -> all genes vs reps (--sensitive, e<=1e-5),
# edges with >=30% coverage of the shorter sequence, connected components, Leiden split of
# components > 7500 genes, reassignment on.
DEFAULT_FAMILY_CONFIG = {"backend": "diamond", "dmnd_cmd": "cluster", "approx_id": 30,
                         "member_cover": 50, "stage2": "cc", "s2_query": "all",
                         "s2_sens": "--sensitive", "s2_evalue": 1e-5, "s2_cov": 0.3,
                         "s2_cov_mode": "shorter", "s2_gmax": 7500, "s2_reassign": True}


# ----------------------------------------------------------------------------
# small utilities
# ----------------------------------------------------------------------------

def md5_file(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _hash(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:16]


def _reset_peak_rss():
    try:
        with open("/proc/self/clear_refs", "w") as fh:
            fh.write("5")
    except OSError:
        pass


def _peak_rss_mb() -> float:
    try:
        with open("/proc/self/status") as fh:
            for line in fh:
                if line.startswith("VmHWM:"):
                    return int(line.split()[1]) / 1024
    except OSError:
        pass
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def _git_commit() -> dict:
    try:
        c = subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True,
                           text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "ancestra"],
                               capture_output=True, text=True).stdout.strip() != ""
        return {"commit": c, "ancestra_dirty": dirty}
    except Exception:  # pragma: no cover
        return {"commit": None}


def _tool_versions() -> dict:
    out = {"python": platform.python_version(), "numpy": np.__version__}
    for tool, args in (("mmseqs", ["version"]), ("diamond", ["version"])):
        try:
            out[tool] = subprocess.run([tool] + args, capture_output=True, text=True).stdout.strip()
        except Exception:
            out[tool] = None
    return out


def atomic_write_text(path: str, text: str) -> None:
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write(text)
    os.replace(tmp, path)


class Runner:
    """Checkpointed stage runner + manifest."""

    def __init__(self, outdir: str, force: bool = False, log=print):
        self.outdir = outdir
        self.sdir = os.path.join(outdir, "stages")
        os.makedirs(self.sdir, exist_ok=True)
        self.force = force
        self.log = log
        self.manifest_path = os.path.join(outdir, "manifest.json")
        self.manifest = {"stages": {}}

    def stage(self, name: str, key: dict, fn):
        d = os.path.join(self.sdir, name)
        h = _hash(key)
        done = os.path.join(d, "DONE")
        if not self.force and os.path.exists(done) and open(done).read().strip() == h:
            with open(os.path.join(d, "stage.json")) as fh:
                info = json.load(fh)
            info["skipped"] = True
            self.manifest["stages"][name] = info
            self.log(f"[ancestra] stage {name}: up to date ({h}), skipped")
            return d, h
        tmp = d + f".tmp{os.getpid()}"
        shutil.rmtree(tmp, ignore_errors=True)
        os.makedirs(tmp)
        _reset_peak_rss()
        t0 = time.time()
        self.log(f"[ancestra] stage {name}: running")
        stats = fn(tmp) or {}
        info = {"key": h, "wall_s": round(time.time() - t0, 3),
                "peak_rss_mb": round(_peak_rss_mb(), 1),
                "children_max_rss_mb": round(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024, 1),
                "stats": stats, "skipped": False}
        with open(os.path.join(tmp, "stage.json"), "w") as fh:
            json.dump(info, fh, indent=1, default=str)
        with open(os.path.join(tmp, "DONE"), "w") as fh:
            fh.write(h + "\n")
        if os.path.exists(d):
            shutil.rmtree(d)
        os.rename(tmp, d)
        self.manifest["stages"][name] = info
        self.log(f"[ancestra] stage {name}: done in {info['wall_s']:.1f}s, "
                 f"peak RSS {info['peak_rss_mb']:.0f} MB")
        return d, h

    def write_manifest(self):
        atomic_write_text(self.manifest_path, json.dumps(self.manifest, indent=1, default=str))


# ----------------------------------------------------------------------------
# families fallback
# ----------------------------------------------------------------------------

def _fallback_families(gt: GeneTable, wd: str, threads: int) -> dict:
    """DIAMOND all-vs-all (sensitive, e<=1e-5, min-cov 50%) + connected components."""
    fa = os.path.join(wd, "genes.fa")
    write_fasta(gt, fa)
    subprocess.run(["diamond", "makedb", "--in", fa, "-d", fa + ".dmnd", "--quiet"], check=True)
    m8 = os.path.join(wd, "ava.m8")
    subprocess.run(["diamond", "blastp", "-q", fa, "-d", fa + ".dmnd", "-o", m8, "--sensitive",
                    "-e", "1e-5", "-k", "500", "--threads", str(threads), "--quiet", "--outfmt", "6",
                    "qseqid", "sseqid", "qcovhsp", "scovhsp"], check=True)
    parent = np.arange(gt.n)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    with open(m8) as fh:
        for line in fh:
            q, t, qc, sc = line.split("\t")
            if min(float(qc), float(sc)) >= 50:
                a, b = find(int(q)), find(int(t))
                if a != b:
                    parent[b] = a
    for f in (fa, fa + ".dmnd", m8):
        os.remove(f)
    return {str(g): str(find(g)) for g in range(gt.n)}


# ----------------------------------------------------------------------------
# pipeline
# ----------------------------------------------------------------------------

def run_pipeline(proteomes: str, tree_path: str, outdir: str, families: str | None = None,
                 fam_cols: tuple = (0, 1), threads: int = 24, scoring: ScoringConfig | None = None,
                 merge: MergeConfig | None = None, family_config: dict | None = None,
                 force: bool = False, log=print) -> dict:
    os.makedirs(outdir, exist_ok=True)
    scoring = scoring or ScoringConfig()
    scoring.threads = threads
    merge = merge or MergeConfig()
    R = Runner(outdir, force=force, log=log)
    t_all = time.time()
    gt = load_proteomes(proteomes)
    tree = load_guide_tree(tree_path, gt.species)
    pfiles = sorted(os.listdir(proteomes))
    inputs = {"proteomes": {f: md5_file(os.path.join(proteomes, f)) for f in pfiles
                            if not f.startswith(".")},
              "tree": md5_file(tree_path),
              "families": md5_file(families) if families else None}
    R.manifest.update(tool="ancestra core v0", argv=sys.argv, git=_git_commit(),
                      versions=_tool_versions(), inputs=inputs,
                      config={"threads": threads, "scoring": asdict(scoring), "merge": asdict(merge),
                              "families": {"file": families, "cols": list(fam_cols)} if families
                              else {"default_clustering": family_config or DEFAULT_FAMILY_CONFIG}},
                      n_genes=gt.n, n_species=len(gt.species))

    # 1. input
    def s_input(d):
        write_gene_table(gt, os.path.join(d, "genes.tsv"))
        atomic_write_text(os.path.join(d, "guide_tree.nwk"), tree.to_newick() + "\n")
        return {"n_genes": gt.n, "n_species": len(gt.species),
                "n_internal_nodes": sum(1 for n in tree.preorder() if n.children)}
    d_in, k_in = R.stage("input", {"inputs": inputs}, s_input)

    # 2. families
    fam_key = {"up": k_in, "file": inputs["families"], "cols": list(fam_cols),
               "cfg": None if families else (family_config or DEFAULT_FAMILY_CONFIG)}

    def s_fam(d):
        info = {}
        if families:
            fam = load_families(families, gt, *fam_cols)
            info["source"] = families
        else:
            try:
                from .families import FamilyConfig, read_families_tsv, run_families
                cfg = FamilyConfig(**(family_config or DEFAULT_FAMILY_CONFIG))
                cfg.threads = threads
                fa = os.path.join(d, "genes.fa")
                write_fasta(gt, fa)
                st = run_families(fa, os.path.join(d, "s0"), cfg)
                m = read_families_tsv(os.path.join(d, "s0", "families.tsv"))
                os.remove(fa)
                info["source"] = "ancestra.families"
                info["s0_timings"] = st.get("timings")
            except ImportError:
                m = _fallback_families(gt, d, threads)
                info["source"] = "fallback diamond+cc"
            fam = families_from_mapping(gt, m)
        np.save(os.path.join(d, "families.npy"), fam)
        sizes = np.bincount(fam)
        info.update(n_families=int(len(sizes)), n_singletons=int((sizes == 1).sum()),
                    max_family=int(sizes.max()), sum_g2=int((sizes.astype(np.int64) ** 2).sum()))
        return info
    d_fam, k_fam = R.stage("families", fam_key, s_fam)
    fam = np.load(os.path.join(d_fam, "families.npy"))

    # 3. scores
    def s_scores(d):
        st, fam2, stats = score_with_rescue(gt, fam, scoring, os.path.join(d, "work"))
        st.save(os.path.join(d, "scores.npz"))
        np.save(os.path.join(d, "families_final.npy"), fam2)
        if st.cross is not None:   # S0-fission signal: hits between different families
            np.savez(os.path.join(d, "cross_family_hits.npz"), **st.cross)
        return stats
    d_sc, k_sc = R.stage("scores", {"up": k_fam, "cfg": scoring.key()}, s_scores)
    st = ScoreTable.load(os.path.join(d_sc, "scores.npz"))

    # 4. merge
    def s_merge(d):
        res = run_merge(gt, tree, st, merge)
        with open(os.path.join(d, "merge.pkl"), "wb") as fh:
            pickle.dump(res, fh, protocol=pickle.HIGHEST_PROTOCOL)
        with open(os.path.join(d, "node_stats.json"), "w") as fh:
            json.dump(res.node_stats, fh, indent=0)
        n_multi = sum(1 for x in res.roots if not isinstance(x, (int, np.integer)))
        return {"n_root_groups": n_multi, "n_root_anc_genes": len(res.roots),
                "merge_s_sum_nodes": round(sum(s["wall_s"] for s in res.node_stats), 3)}
    d_m, k_m = R.stage("merge", {"up": k_sc, "cfg": merge.key()}, s_merge)
    with open(os.path.join(d_m, "merge.pkl"), "rb") as fh:
        res: MergeResult = pickle.load(fh)
    R.manifest["nodes"] = [{k: s[k] for k in ("node", "n_children", "n_units", "n_anc_out",
                                              "n_merged_groups", "n_edges", "wall_s")}
                           for s in res.node_stats]

    # 5. output
    def s_out(d):
        ox = to_orthoxml(res, gt, tree)
        x = os.path.join(d, "hogs.orthoxml")
        write_orthoxml(ox, x)
        ok, errs = validate_xsd(x)
        if not ok:
            raise RuntimeError("OrthoXML failed XSD validation: " + "; ".join(errs[:5]))
        pairs = hogs_to_pairs(ox, with_level=True)
        with open(os.path.join(d, "pairs.tsv"), "w") as fh:
            fh.write("#gene_a\tgene_b\tlevel\n")
            for (a, b), lvl in sorted(pairs.items()):
                fh.write(f"{a}\t{b}\t{lvl}\n")
        ogs = to_orthobench_ogs(ox, None, os.path.join(d, "orthobench_ogs.txt"))
        return {"xsd_valid": ok, "n_top_groups": len(ox.groups), "n_pairs": len(pairs),
                "n_ogs_multi": sum(1 for o in ogs if len(o) > 1)}
    d_o, k_o = R.stage("output", {"up": k_m}, s_out)
    for f in ("hogs.orthoxml", "pairs.tsv", "orthobench_ogs.txt"):
        src, dst = os.path.join(d_o, f), os.path.join(outdir, f)
        tmp = dst + ".tmp"
        shutil.copyfile(src, tmp)
        os.replace(tmp, dst)
    R.manifest["total_wall_s"] = round(time.time() - t_all, 3)
    R.manifest["outputs"] = {f: md5_file(os.path.join(outdir, f))
                             for f in ("hogs.orthoxml", "pairs.tsv", "orthobench_ogs.txt")}
    R.write_manifest()
    return R.manifest


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------

def build_parser(ap: argparse.ArgumentParser | None = None) -> argparse.ArgumentParser:
    ap = ap or argparse.ArgumentParser(prog="python -m ancestra run", description=__doc__.split("\n")[0])
    ap.add_argument("--proteomes", required=True, help="directory with one FASTA per species")
    ap.add_argument("--tree", required=True, help="rooted Newick guide tree (internal names = levels)")
    ap.add_argument("--families", help="gene -> family TSV (default: S0 clustering)")
    ap.add_argument("--families-cols", default="0,1", help="gene,family column indices (0-based)")
    ap.add_argument("--family-config", help="JSON with ancestra.families.FamilyConfig fields")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=24)
    ap.add_argument("--backend", default="hybrid", choices=["hybrid", "diamond", "mmseqs"])
    ap.add_argument("--arm", default=None, choices=sorted(ARMS),
                    help="named merge configuration (default: MergeConfig defaults)")
    ap.add_argument("--config", help="JSON {'scoring': {...}, 'merge': {...}} overriding defaults")
    ap.add_argument("--force", action="store_true", help="recompute all stages")
    return ap


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    sc = ScoringConfig(backend=a.backend)
    mc = MergeConfig(**asdict(ARMS[a.arm])) if a.arm else MergeConfig()
    if a.config:
        with open(a.config) as fh:
            c = json.load(fh)
        for k, v in c.get("scoring", {}).items():
            setattr(sc, k, v)
        for k, v in c.get("merge", {}).items():
            setattr(mc, k, v)
    fc = None
    if a.family_config:
        with open(a.family_config) as fh:
            fc = json.load(fh)
    cols = tuple(int(x) for x in a.families_cols.split(","))
    m = run_pipeline(a.proteomes, a.tree, a.out, families=a.families, fam_cols=cols,
                     threads=a.threads, scoring=sc, merge=mc, family_config=fc, force=a.force,
                     log=lambda s: print(s, file=sys.stderr))
    print(f"[ancestra] done in {m['total_wall_s']:.1f}s -> {a.out}/hogs.orthoxml", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
