#!/usr/bin/env python
"""S8 driver: apparent-loss candidates, controls, genome-wide search, known biology.

Per target species (``--species``), writes into ``--out``:
  cand_<sp>.tsv      locus classification of every apparent-loss candidate
  pos_<sp>.tsv       positive controls (complete single-copy HOGs, target gene hidden)
  mask_<sp>.tsv      the positive-control loci with the hidden gene span set to N
  neg_<sp>.tsv       random windows (size of a candidate locus) with candidate queries
  gw_<sp>.tsv        genome-wide miniprot for B/C/D candidates (relocated copies)
  known_<sp>.tsv     GULO / UOX direct tests (human, chimp)
Summaries: scripts/summarize_S8.py.
"""
import argparse
import json
import os
import pickle
import random
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ancestra.geneloss import (MAMMALS, CdsIndex, Locus, SOURCE_ORDER, FastaIndex, GeneOrderIndex, HogIndex,  # noqa: E402
                               LocusSearcher, parse_miniprot_paf, random_window, read_fasta,
                               select_apparent_losses, select_complete_single_copy, syntenic_locus,
                               write_fasta)

GENOME = {
    "Homo_sapiens": "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
    "Mus_musculus": "Mus_musculus.GRCm38.dna_sm.primary_assembly.fa",
    "Rattus_norvegicus": "Rattus_norvegicus.Rnor_6.0.dna_sm.toplevel.fa",
    "Pan_troglodytes": "Pan_troglodytes.Pan_tro_3.0.dna_sm.toplevel.fa",
    "Canis_familiaris": "Canis_familiaris.CanFam3.1.dna_sm.toplevel.fa",
    "Monodelphis_domestica": "Monodelphis_domestica.ASM229v1.dna_sm.toplevel.fa",
}
PROTEOME = {
    "Homo_sapiens": "Homo_sapiens.GRCh38.pep.all.fa", "Mus_musculus": "Mus_musculus.GRCm38.pep.all.fa",
    "Rattus_norvegicus": "Rattus_norvegicus.Rnor_6.0.pep.all.fa", "Pan_troglodytes": "Pan_troglodytes.Pan_tro_3.0.pep.all.fa",
    "Canis_familiaris": "Canis_familiaris.CanFam3.1.pep.all.fa",
    "Monodelphis_domestica": "Monodelphis_domestica.ASM229v1.pep.all.fa",
}
OUTGROUP = ["Gallus_gallus", "Danio_rerio", "Tetraodon_nigroviridis"]
KNOWN = {  # name -> mouse protein (query HOG), expected targets
    "GULO": ("ENSMUSP00000060912", ["Homo_sapiens", "Pan_troglodytes"]),
    "UOX": ("ENSMUSP00000029837", ["Homo_sapiens", "Pan_troglodytes"]),
}


def load_hi(path_xml, tree, cache):
    if os.path.exists(cache):
        return pickle.load(open(cache, "rb"))
    hi = HogIndex.from_orthoxml(path_xml, tree, "Theria", "Amniota")
    pickle.dump(hi, open(cache, "wb"))
    return hi


def gff_overlap(gff, contig, blocks):
    """Ensembl genes (any biotype) overlapping the hit's CDS blocks, by overlap bp."""
    if gff is None or not blocks:
        return "", ""
    g = gff[gff.contig == contig]
    lo, hi = min(x for x, _ in blocks), max(y for _, y in blocks)
    g = g[(g.start <= hi) & (g.end >= lo)]
    if len(g) == 0:
        return "", ""
    ov = [sum(max(0, min(e, ge) - max(s, gs) + 1) for s, e in blocks) for gs, ge in zip(g.start, g.end)]
    g = g.assign(ov=ov)
    g = g[g.ov > 0].sort_values("ov", ascending=False)
    return ";".join(g.name.replace("", "-").astype(str).head(3)), ";".join(g.biotype.head(3))


def parse_blocks(txt):
    if not isinstance(txt, str) or not txt:
        return []
    return [tuple(int(v) for v in x.split("-")) for x in txt.split(",")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", nargs="+", default=MAMMALS)
    ap.add_argument("--orthoxml", default="work/core_v0_ob/run_v2/hogs.orthoxml")
    ap.add_argument("--tree", default="results/S1_speciestree/ob/ncbi.nwk")
    ap.add_argument("--tables", default="data/genomes/orthobench")
    ap.add_argument("--proteomes", default="data/orthobench/Open_Orthobench/BENCHMARKS/Input")
    ap.add_argument("--genomes", default="data/genomes_dna")
    ap.add_argument("--work", default="work/S8")
    ap.add_argument("--out", default="results/S8_geneloss")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--n-pos", type=int, default=60)
    ap.add_argument("--n-neg", type=int, default=100)
    ap.add_argument("--max-cand", type=int, default=0, help="0 = all")
    ap.add_argument("--seed", type=int, default=8)
    ap.add_argument("--pad", type=int, default=200_000)
    ap.add_argument("--skip-gw", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    os.makedirs(a.work, exist_ok=True)
    t0 = time.time()
    hi = load_hi(a.orthoxml, a.tree, os.path.join(a.work, "hogindex.pkl"))
    go = GeneOrderIndex.from_dir(a.tables, MAMMALS)
    prot = {}
    for sp in MAMMALS:
        prot.update(read_fasta(os.path.join(a.proteomes, PROTEOME[sp])))
    n_exons = {}
    for sp in MAMMALS:
        n_exons.update(pd.read_parquet(os.path.join(a.work, f"cds_{sp}.parquet")).groupby("protein_id").size().to_dict())
    cand = select_apparent_losses(hi, MAMMALS, max_copies=2, outgroup=OUTGROUP,
                                  root_species={"Monodelphis_domestica"})
    cand.to_csv(os.path.join(a.out, "candidates_all.tsv"), sep="\t", index=False)
    pos_hogs = select_complete_single_copy(hi, MAMMALS)
    print(f"[S8] loaded in {time.time()-t0:.0f}s; candidates {len(cand)}; complete 1:1 HOGs {len(pos_hogs)}", flush=True)

    for T in a.species:
        ts = time.time()
        gpath = os.path.join(a.genomes, GENOME[T])
        if not os.path.exists(gpath + ".fai"):
            print(f"[S8] {T}: genome not ready, skip", flush=True)
            continue
        genome = FastaIndex(gpath)
        clen = genome.lengths()
        gffp = os.path.join(a.work, f"gff_genes_{T}.parquet")
        gff = pd.read_parquet(gffp) if os.path.exists(gffp) else None
        cds = CdsIndex.from_parquet(os.path.join(a.work, f"cds_{T}.parquet"))
        S = LocusSearcher(hi, go, T, genome, prot, os.path.join(a.work, "tmp", T), pad=a.pad, cds=cds,
                          n_exons=n_exons)
        rng = random.Random(f"{a.seed}-{T}")
        srcs = SOURCE_ORDER[T]

        def annotate(rec):
            if rec.get("hit_start") is not None and rec.get("contig"):
                rec["gff_names"], rec["gff_biotypes"] = gff_overlap(gff, rec["contig"], parse_blocks(rec.get("hit_blocks")))
            return rec

        # --- candidates
        cT = cand[cand.target == T]
        if a.max_cand:
            cT = cT.head(a.max_cand)
        jobs = []
        for r in cT.itertuples():
            loc = syntenic_locus(hi, go, r.hog, T, srcs, pad=a.pad, contig_len=clen)
            jobs.append((f"cand{r.hog}", r.hog, loc, frozenset(), None, None))

        # --- positive controls (hidden target gene) + masked
        pos = rng.sample(pos_hogs, min(a.n_pos, len(pos_hogs)))
        pjobs, mjobs = [], []
        for h in pos:
            tg = hi.in_species(h, T)[0]
            loc = syntenic_locus(hi, go, h, T, srcs, pad=a.pad, contig_len=clen)
            pjobs.append((f"pos{h}", h, loc, frozenset([tg]), None, tg))
            r = go.row(T, tg)
            if r is not None and loc is not None:
                mjobs.append((f"mask{h}", h, loc, frozenset([tg]), (int(r["start"]), int(r["end"])), None))

        # --- negative controls: random windows, candidate queries
        njobs = []
        allowed = set(go.by_contig[T].keys())
        sized = [j for j in jobs if j[2] is not None]
        for j in rng.sample(sized, min(a.n_neg, len(sized))):
            loc = j[2]
            c, s, e = random_window(genome, loc.end - loc.start + 1, rng, allowed=allowed)
            nl = Locus(T, c, s, e, "bracketed", loc.source, loc.source_gene)
            njobs.append((f"neg{j[1]}", j[1], nl, frozenset(), None, None))

        def run(job):
            jid, h, loc, hide, mask, keep = job
            try:
                return annotate(S.search(jid, h, loc, hide=hide, mask=mask, keep_hidden_as_fission=keep))
            except Exception as ex:  # keep going, record the failure
                return {"job": jid, "hog": h, "target": T, "cls": "ERROR", "error": repr(ex)[:200]}

        for name, jl in (("cand", jobs), ("pos", pjobs), ("mask", mjobs), ("neg", njobs)):
            with ThreadPoolExecutor(a.threads) as ex:
                recs = list(ex.map(run, jl))
            df = pd.DataFrame(recs)
            if name == "pos":
                tgt = {f"pos{h}": hi.in_species(h, T)[0] for h in pos}
                df["hidden_gene"] = df.job.map(tgt)
                df["recovered"] = (df.cls == "A_annotated") & (df.overlap_gene == df.hidden_gene)
                # as if the hidden gene were unannotated: would the A_gap rule fire?
                df["recovered_gap_rule"] = df.recovered & (df.nq_intact > df.nq_disrupted)
            df.to_csv(os.path.join(a.out, f"{name}_{T}.tsv"), sep="\t", index=False)
            print(f"[S8] {T} {name}: n={len(df)} {df.cls.value_counts().to_dict()} ({time.time()-ts:.0f}s)", flush=True)

        # --- known biology (GULO, UOX) with the HOG of the mouse protein
        krecs = []
        for nm, (mp, targets) in KNOWN.items():
            if T not in targets or mp not in hi.hog_of:
                continue
            h = hi.hog_of[mp]
            loc = syntenic_locus(hi, go, h, T, srcs, pad=a.pad, contig_len=clen)
            hide = frozenset(hi.in_species(h, T))
            rec = annotate(S.search(f"known_{nm}", h, loc, hide=hide))
            rec["gene"] = nm
            rec["hog_species"] = ",".join(sorted({hi.species_of[g] for g in hi.members[h]}))
            rec["in_candidates"] = bool(((cand.hog == h)).any())
            krecs.append(rec)
        if krecs:
            pd.DataFrame(krecs).to_csv(os.path.join(a.out, f"known_{T}.tsv"), sep="\t", index=False)
            print(f"[S8] {T} known: {[(r['gene'], r['cls'], r.get('fs'), r.get('st'), r.get('qcov')) for r in krecs]}", flush=True)

        # --- genome-wide miniprot for B/C/D candidates
        if a.skip_gw:
            continue
        cd = pd.read_csv(os.path.join(a.out, f"cand_{T}.tsv"), sep="\t")
        bcd = cd[cd.cls.isin(["B_pseudo", "B_partial", "C_deleted", "D_undetermined"])]
        q = {}
        for h in bcd.hog:
            for g, s in S.queries(int(h)).items():
                q[g] = s
        if not q:
            continue
        idx = os.path.join(a.work, f"{T}.mpi")
        if not os.path.exists(idx):
            subprocess.run(["miniprot", "-t", str(a.threads), "-d", idx, gpath], check=True,
                           stderr=subprocess.DEVNULL)
        qf = os.path.join(a.work, f"gw_queries_{T}.faa")
        write_fasta(qf, q)
        r = subprocess.run(["miniprot", "--gff", "-t", str(a.threads), "--outs=0.5", "-p", "0.5", "-N", "20", idx, qf],
                           capture_output=True, text=True, check=True)
        open(os.path.join(a.work, f"gw_{T}.paf"), "w").write(r.stdout)
        hits = parse_miniprot_paf(r.stdout)
        byq = {}
        for hh in hits:
            byq.setdefault(hh.query, []).append(hh)
        grec = []
        for row in bcd.itertuples():
            h = int(row.hog)
            best, best_any, npar = None, None, 0
            for g in S.queries(h):
                for hh in byq.get(g, []):
                    fis, par = S.paralog_or_fission(hh.contig, hh.blocks, h)
                    if par is not None and fis is None:
                        npar += 1
                        continue
                    if best_any is None or (hh.qcov, hh.score) > (best_any[0].qcov, best_any[0].score):
                        best_any = (hh, fis)
                    if hh.intact and (best is None or hh.score > best[0].score):
                        best = (hh, fis)
            rec = {"hog": h, "target": T, "locus_cls": row.cls, "gw_paralog_hits": npar}
            for tag, b in (("intact", best), ("any", best_any)):
                if b is None:
                    rec[f"gw_{tag}"] = False
                    continue
                hh, fis = b
                inside = (isinstance(row.contig, str) and hh.contig == row.contig
                          and hh.ts + 1 <= row.end and hh.te >= row.start)
                gn, gb = gff_overlap(gff, hh.contig, hh.blocks)
                rec.update({f"gw_{tag}": True, f"gw_{tag}_contig": hh.contig, f"gw_{tag}_start": hh.ts + 1,
                            f"gw_{tag}_end": hh.te, f"gw_{tag}_qcov": round(hh.qcov, 3),
                            f"gw_{tag}_identity": round(hh.identity, 3), f"gw_{tag}_fs": hh.fs, f"gw_{tag}_st": hh.st,
                            f"gw_{tag}_n_introns": hh.n_introns, f"gw_{tag}_query": hh.query,
                            f"gw_{tag}_in_locus": bool(inside), f"gw_{tag}_annotated": fis or "",
                            f"gw_{tag}_gff": gn, f"gw_{tag}_biotype": gb})
            grec.append(rec)
        pd.DataFrame(grec).to_csv(os.path.join(a.out, f"gw_{T}.tsv"), sep="\t", index=False)
        g = pd.DataFrame(grec)
        print(f"[S8] {T} genome-wide: n={len(g)} intact_elsewhere={int(g.gw_intact.sum())} ({time.time()-ts:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
