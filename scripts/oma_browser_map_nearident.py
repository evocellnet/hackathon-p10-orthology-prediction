#!/usr/bin/env python3
"""Route 6 of the OrthoBench -> OMA mapping: near-identical sequence match for proteins WITHOUT an in-HOG OMA entry
after routes 1-5 (mostly chicken: Ensembl 99 GRCg6a vs OMA's RefSeq GRCg7b). DIAMOND blastp (same species only) against the
OMA entries of that species that are HOG members and not already used by another OrthoBench protein; accept the best hit
iff pident >= 95, query and subject coverage >= 90 %, unique best bitscore; an OMA entry claimed by several queries goes to the
highest bitscore. Usage: oma_browser_map_nearident.py mapping.tsv hog_paths.tsv oma-seqs.fa.gz <OB Input dir> workdir out.tsv threads
"""
import sys, os, gzip, glob, subprocess, collections
mp, paths, seqs, obdir, wd, out, thr = sys.argv[1:8]
os.makedirs(wd, exist_ok=True)
rows = [l.rstrip("\n").split("\t") for l in open(mp)]
hdr, rows = rows[0], rows[1:]
SPC = {"Caenorhabditis_elegans": "CAEEL", "Canis_familiaris": "CANLF", "Ciona_intestinalis": "CIOIN", "Danio_rerio": "DANRE",
       "Drosophila_melanogaster": "DROME", "Gallus_gallus": "CHICK", "Homo_sapiens": "HUMAN", "Monodelphis_domestica": "MONDO",
       "Mus_musculus": "MOUSE", "Pan_troglodytes": "PANTR", "Rattus_norvegicus": "RATNO", "Tetraodon_nigroviridis": "TETNG"}
in_hog = set(l.split("\t", 1)[0] for l in open(paths))
used = set(r[2] for r in rows if r[4] == "inhog")
todo = [r for r in rows if r[4] != "inhog"]
need_sp = set(SPC[r[1]] for r in todo)
def read_fa(fh):
    h = None; s = []
    for l in fh:
        if l.startswith(">"):
            if h is not None: yield h, "".join(s)
            h = l[1:].split()[0]; s = []
        else: s.append(l.strip())
    if h is not None: yield h, "".join(s)
ob = {}
for f in glob.glob(os.path.join(obdir, "*.fa")):
    for h, s in read_fa(open(f)): ob[h] = s.rstrip("*")
with open(f"{wd}/db.fa", "w") as fo:
    for h, s in read_fa(gzip.open(seqs, "rt")):
        if h[:5] in need_sp and h in in_hog and h not in used: fo.write(f">{h}\n{s}\n")
with open(f"{wd}/q.fa", "w") as fo:
    for r in todo: fo.write(f">{r[0]}|{SPC[r[1]]}\n{ob[r[0]]}\n")
subprocess.run(["diamond", "makedb", "--in", f"{wd}/db.fa", "-d", f"{wd}/db", "--quiet"], check=True)
subprocess.run(["diamond", "blastp", "-q", f"{wd}/q.fa", "-d", f"{wd}/db", "-o", f"{wd}/hits.tsv", "--threads", thr, "--quiet",
                "--id", "95", "--query-cover", "90", "--subject-cover", "90", "-k", "5", "--ultra-sensitive",
                "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", "scovhsp", "bitscore"], check=True)
best = collections.defaultdict(list)
for l in open(f"{wd}/hits.tsv"):
    q, s, pid, qc, sc, bs = l.split("\t"); qid, code = q.split("|")
    if s[:5] != code: continue
    best[qid].append((float(bs), s))
claim = {}
for q, hs in best.items():
    hs.sort(reverse=True)
    if len(hs) > 1 and hs[0][0] == hs[1][0]: continue  # ambiguous
    bs, s = hs[0]
    if s not in claim or claim[s][0] < bs: claim[s] = (bs, q)
q2s = {q: s for s, (bs, q) in claim.items()}
n = collections.Counter()
with open(out, "w") as fo:
    fo.write("\t".join(hdr) + "\n")
    for r in rows:
        if r[4] != "inhog" and r[0] in q2s:
            r = [r[0], r[1], q2s[r[0]], "near_identical_diamond", "inhog"]; n[r[1]] += 1
        fo.write("\t".join(r) + "\n")
print("near-identical rescues:", sum(n.values()), dict(n))
