#!/usr/bin/env python3
"""List (and optionally download) the public ortholog-prediction files submitted by the
participants of the QfO 2022 benchmark (OpenEBench event OEBE0020000003). The files are
deposited on EUDAT B2SHARE (CC-BY 4.0) with titles
"Quest for Orthologs - 2022 Reference Proteomes participant dataset submitted by <tool>".
(The DOIs stored in OpenEBench point to legacy B2SHARE ids that now 404, so we search by title.)

Usage: fetch_qfo_participant_predictions.py [--list out.tsv] [--download TOOL ...] [--outdir DIR]
  TOOL is matched exactly against the participant name, e.g. FastOMA_v0.1.4_nov OrthoFinder
"""
import argparse, hashlib, json, os, sys, urllib.parse, urllib.request

TITLE = "Quest for Orthologs - 2022 Reference Proteomes participant dataset submitted by "
p = argparse.ArgumentParser()
p.add_argument("--list", default=None)
p.add_argument("--download", nargs="*", default=[])
p.add_argument("--outdir", default=os.path.join(os.environ.get("DATA_DIR", "data"), "qfo", "participants_2022"))
a = p.parse_args()
url = "https://b2share.eudat.eu/api/records?size=100&q=" + urllib.parse.quote('"%s"' % TITLE.strip())
hits = json.load(urllib.request.urlopen(url, timeout=300))["hits"]["hits"]
rows = []
for h in hits:
    title = h["metadata"].get("title", "")
    if not title.startswith(TITLE):
        continue
    tool = title[len(TITLE):]
    rec = json.load(urllib.request.urlopen(f"https://b2share.eudat.eu/api/records/{h['id']}", timeout=300))
    for k, f in rec.get("files", {}).get("entries", {}).items():
        rows.append((tool, h["id"], k, f["size"], f["checksum"], f["links"]["content"], h["metadata"].get("publication_date", "")))
rows.sort()
out = open(a.list, "w") if a.list else sys.stdout
print("tool\trecord\tfile\tsize\tchecksum\turl\tpublication_date", file=out)
for r in rows:
    print("\t".join(map(str, r)), file=out)
for tool in a.download:
    os.makedirs(a.outdir, exist_ok=True)
    sel = [r for r in rows if r[0] == tool]
    if not sel:
        print("no record for", tool, file=sys.stderr)
    for r in sel:
        fn = os.path.join(a.outdir, r[2])
        if not os.path.exists(fn):
            print("downloading", r[2], r[3], "bytes", file=sys.stderr)
            with urllib.request.urlopen(r[5], timeout=600) as src, open(fn + ".part", "wb") as dst:
                while True:
                    b = src.read(1 << 22)
                    if not b: break
                    dst.write(b)
            os.rename(fn + ".part", fn)
        h = hashlib.md5()
        with open(fn, "rb") as f:
            for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
        print(fn, "md5", h.hexdigest(), "OK" if r[4].endswith(h.hexdigest()) else "MISMATCH " + r[4], file=sys.stderr)
