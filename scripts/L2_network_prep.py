#!/usr/bin/env python
"""L2-network data preparation: STRING v12 own-species channels and PaxDb v5 abundance,
mapped onto the 12 OrthoBench proteomes.

Circularity rules (docs/L2_network.md): only the organism's OWN channels are kept
(``experiments``, ``coexpression``, ``database``; never ``*_transferred``, ``textmining``,
``homology``, ``neighborhood``, ``fusion`` or ``cooccurence``, which are either transferred by
orthology or genomic-context evidence). ``database`` is kept only as a flagged secondary
channel: KEGG/Reactome pathway membership in non-human species is itself assigned by orthology.

ID mapping STRING protein -> OrthoBench protein, in priority order:
  1. Ensembl protein ID (STRING ID without taxid prefix, or an Ensembl_protein_id /
     Ensembl_translation alias) equal to the OrthoBench ID (version-stripped);
  2. identical sequence (md5);
  3. same gene (STRING ``Ensembl_gene`` / ``WormBase_gene`` alias == OrthoBench gene_id from the
     gene-order tables, data/genomes/orthobench).
Output (MAIN/work/l2_network/): ``edges_<species>.parquet`` (a, b = OrthoBench IDs, a<b,
exp, coexp, db channel scores, max over STRING proteins mapping to the pair) and
``abund_<species>.parquet`` (OrthoBench ID, ppm); ``mapping_coverage.json``.
"""
import gzip, hashlib, json, os, sys, time

import numpy as np
import pandas as pd

MAIN = os.environ.get("P10_MAIN", "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction")
OB_IN = f"{MAIN}/data/orthobench/Open_Orthobench/BENCHMARKS/Input"
SD = f"{MAIN}/data/extra/l2_network/string"
PD = f"{MAIN}/data/extra/l2_network/paxdb"
W = f"{MAIN}/work/l2_network"

SPECIES = {  # OrthoBench species -> (fasta, NCBI taxid)
    "Caenorhabditis_elegans": ("Caenorhabditis_elegans.WBcel235.pep.all.fa", 6239),
    "Canis_familiaris": ("Canis_familiaris.CanFam3.1.pep.all.fa", 9615),
    "Ciona_intestinalis": ("Ciona_intestinalis.KH.pep.all.fa", 7719),
    "Danio_rerio": ("Danio_rerio.GRCz11.pep.all.fa", 7955),
    "Drosophila_melanogaster": ("Drosophila_melanogaster.BDGP6.28.pep.all.fa", 7227),
    "Gallus_gallus": ("Gallus_gallus.GRCg6a.pep.all.fa", 9031),
    "Homo_sapiens": ("Homo_sapiens.GRCh38.pep.all.fa", 9606),
    "Monodelphis_domestica": ("Monodelphis_domestica.ASM229v1.pep.all.fa", 13616),
    "Mus_musculus": ("Mus_musculus.GRCm38.pep.all.fa", 10090),
    "Pan_troglodytes": ("Pan_troglodytes.Pan_tro_3.0.pep.all.fa", 9598),
    "Rattus_norvegicus": ("Rattus_norvegicus.Rnor_6.0.pep.all.fa", 10116),
    "Tetraodon_nigroviridis": ("Tetraodon_nigroviridis.TETRAODON8.pep.all.fa", 99883),
}
PAXDB = {6239: "6239-WHOLE_ORGANISM-integrated.txt", 9615: "9615-GPM_2012_09_Canis_familiaris.txt",
         7955: "7955-WHOLE_ORGANISM-integrated.txt", 7227: "7227-WHOLE_ORGANISM-integrated.txt",
         9031: "9031-WHOLE_ORGANISM-integrated.txt", 9606: "9606-WHOLE_ORGANISM-integrated.txt",
         10090: "10090-WHOLE_ORGANISM-integrated.txt", 9598: "9598-Chimp_iBAQ_Khan_2013.txt",
         10116: "10116-WHOLE_ORGANISM-integrated.txt"}
OWN = ["experiments", "coexpression", "database"]


def read_fasta(path, opener=open):
    seqs, name, buf = {}, None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def md5(s):
    return hashlib.md5(s.rstrip("*").upper().encode()).hexdigest()


def strip_v(x):
    return x.rsplit(".", 1)[0] if x.count(".") == 1 and x.rsplit(".", 1)[1].isdigit() else x


def build_map(sp, fa, taxid):
    ob = read_fasta(f"{OB_IN}/{fa}")
    ob_ids = list(ob)
    ob_by_id = {strip_v(i): i for i in ob_ids}
    ob_by_md5 = {}
    for i, s in ob.items():
        ob_by_md5.setdefault(md5(s), []).append(i)
    go = pd.read_parquet(f"{MAIN}/data/genomes/orthobench/{sp}.parquet", columns=["protein_id", "gene_id"])
    ob_by_gene = {}
    for p, g in zip(go.protein_id, go.gene_id):
        ob_by_gene.setdefault(strip_v(str(g)), []).append(p)
    st = read_fasta(f"{SD}/{taxid}.protein.sequences.v12.0.fa.gz", gzip.open)
    al = pd.read_csv(f"{SD}/{taxid}.protein.aliases.v12.0.txt.gz", sep="\t", comment=None, header=0,
                     names=["sid", "alias", "source"], dtype=str)
    al_prot = al[al.source.isin(["Ensembl_protein_id", "Ensembl_translation"])]
    al_gene = al[al.source.isin(["Ensembl_gene", "WormBase_gene", "Ensembl_wormbase_gene"])]
    prot_alias = {}
    for s_, a in zip(al_prot.sid, al_prot.alias):
        prot_alias.setdefault(s_, set()).add(strip_v(a))
    gene_alias = {}
    for s_, a in zip(al_gene.sid, al_gene.alias):
        gene_alias.setdefault(s_, set()).add(strip_v(a))
    m, how = {}, {"id": 0, "md5": 0, "gene": 0, "none": 0}
    for sid, seq in st.items():
        core = sid.split(".", 1)[1]
        cands = [core, strip_v(core)] + sorted(prot_alias.get(sid, ()))
        hit = next((ob_by_id[c] for c in cands if c in ob_by_id), None)
        if hit is not None:
            m[sid] = [hit]; how["id"] += 1; continue
        h = ob_by_md5.get(md5(seq))
        if h:
            m[sid] = h; how["md5"] += 1; continue
        gg = [p for g in sorted(gene_alias.get(sid, ())) for p in ob_by_gene.get(g, [])]
        if gg:
            m[sid] = sorted(set(gg)); how["gene"] += 1; continue
        how["none"] += 1
    mapped_ob = {p for v in m.values() for p in v}
    cov = dict(n_ob=len(ob_ids), n_string=len(st), string_mapped_by=how,
               ob_mapped=len(mapped_ob), ob_mapped_frac=round(len(mapped_ob) / len(ob_ids), 4))
    # alias -> STRING ids (for PaxDb files keyed by older STRING IDs, e.g. worm sequence names)
    a2s = {}
    for s_, a in zip(al.sid, al.alias):
        a2s.setdefault(a, set()).add(s_)
    return m, cov, a2s


def edges(sp, taxid, m):
    t0 = time.time()
    d = pd.read_csv(f"{SD}/{taxid}.protein.links.full.v12.0.txt.gz", sep=" ",
                    usecols=["protein1", "protein2"] + OWN)
    d = d[(d[OWN] > 0).any(axis=1)]
    d = d[d.protein1 < d.protein2]          # each undirected link once (file lists both)
    rows = []
    for col_a, col_b in (("protein1", "protein2"),):
        a = d[col_a].map(m); b = d[col_b].map(m)
        ok = a.notna() & b.notna()
        x = d[ok].assign(a=a[ok], b=b[ok]).explode("a").explode("b")
        rows.append(x)
    x = pd.concat(rows)
    x = x[x.a != x.b]
    lo = np.where(x.a < x.b, x.a, x.b); hi = np.where(x.a < x.b, x.b, x.a)
    x = pd.DataFrame({"a": lo, "b": hi, "exp": x.experiments.values, "coexp": x.coexpression.values,
                      "db": x.database.values})
    x = x.groupby(["a", "b"], as_index=False).max()
    x.to_parquet(f"{W}/edges_{sp}.parquet", index=False)
    st = {f"n_edges_{c}_ge{t}": int((x[c] >= t).sum()) for c in ("exp", "coexp", "db") for t in (150, 400, 700)}
    for c in ("exp", "coexp", "db"):
        y = x[x[c] >= 400]
        st[f"n_nodes_{c}_ge400"] = int(len(set(y.a) | set(y.b)))
    y = x[(x.exp >= 400) | (x.coexp >= 400)]
    st["n_nodes_expcoexp_ge400"] = int(len(set(y.a) | set(y.b)))
    st["n_edges_expcoexp_ge400"] = int(len(y))
    st["wall_s"] = round(time.time() - t0, 1)
    return st


def abundance(sp, taxid, m, a2s):
    if taxid not in PAXDB:
        return None
    d = pd.read_csv(f"{PD}/{PAXDB[taxid]}", sep="\t", comment="#", header=None, names=["sid", "ppm"],
                    usecols=[0, 1])
    d0 = d.copy()
    def look(sid):
        if sid in m:
            return m[sid]
        core = sid.split(".", 1)[1]
        hits = sorted({p for s_ in a2s.get(core, ()) for p in m.get(s_, [])})
        return hits or None
    d["ob"] = d.sid.map(look)
    d = d[d.ob.notna()].explode("ob")
    d = d.groupby("ob", as_index=False).ppm.max()
    d.to_parquet(f"{W}/abund_{sp}.parquet", index=False)
    return dict(paxdb_file=PAXDB[taxid], n_paxdb=int(len(d0)), n_ob_with_abundance=int(len(d)), n_pos=int((d.ppm > 0).sum()))


def humap():
    """hu.MAP 3.0 (human CF-MS + AP-MS integration): complex co-membership pairs, UniProt AC ->
    STRING (UniProt_AC alias) -> OrthoBench. Written to humap_Homo_sapiens.parquet (a, b, prob)."""
    fa, taxid = SPECIES["Homo_sapiens"]
    m, c, a2s = build_map("Homo_sapiens", fa, taxid)
    d = pd.read_csv(f"{MAIN}/data/extra/l2_network/humap3/humap3_complex_pairs_20240929.pairsWprob", sep="\t",
                    header=None, names=["u1", "u2", "prob"])
    look = lambda u: sorted({p for s_ in a2s.get(u, ()) for p in m.get(s_, [])}) or None
    d["a"] = d.u1.map(look); d["b"] = d.u2.map(look)
    n0 = len(d); d = d[d.a.notna() & d.b.notna()].explode("a").explode("b")
    d = d[d.a != d.b]
    lo = np.where(d.a < d.b, d.a, d.b); hi = np.where(d.a < d.b, d.b, d.a)
    x = pd.DataFrame({"a": lo, "b": hi, "prob": d.prob.values}).groupby(["a", "b"], as_index=False).max()
    x.to_parquet(f"{W}/humap_Homo_sapiens.parquet", index=False)
    st = dict(n_pairs_file=n0, n_pairs_mapped=int(len(x)), n_nodes=int(len(set(x.a) | set(x.b))))
    print("humap", st)
    return st


def main():
    os.makedirs(W, exist_ok=True)
    if sys.argv[1:] == ["--humap"]:
        covf = f"{W}/mapping_coverage.json"
        cov = json.load(open(covf)); cov["humap3_Homo_sapiens"] = humap()
        json.dump(cov, open(covf, "w"), indent=1)
        return
    only = sys.argv[1:] or list(SPECIES)
    covf = f"{W}/mapping_coverage.json"
    cov = json.load(open(covf)) if os.path.exists(covf) else {}
    for sp in only:
        fa, taxid = SPECIES[sp]
        t0 = time.time()
        m, c, a2s = build_map(sp, fa, taxid)
        c["edges"] = edges(sp, taxid, m)
        c["abundance"] = abundance(sp, taxid, m, a2s)
        c["wall_s"] = round(time.time() - t0, 1)
        cov[sp] = c
        print(sp, json.dumps(c), flush=True)
        json.dump(cov, open(covf, "w"), indent=1)


if __name__ == "__main__":
    main()
