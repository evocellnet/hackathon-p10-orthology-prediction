#!/usr/bin/env python3
"""Harvest public QfO benchmark (2022 reference proteomes) assessment metrics from OpenEBench.

Source: OpenEBench REST API, community OEBC002 (Quest for Orthologs).
  https://openebench.bsc.es/api/scientific/staged/Dataset            (list of dataset ids)
  https://openebench.bsc.es/api/scientific/staged/Dataset/<id>       (one dataset)
Assessment datasets carry inline {"value", "error"}; their `name` encodes metric/challenge/participant.
The 2022 event challenges are OEBX0020000015 .. OEBX002000001H.

Usage: python fetch_qfo_public_results.py OUT.tsv   (needs internet; ~4.5k HTTP requests, ~2-5 min)
Output: long-format TSV participant, challenge, metric, value, error, dataset_id
"""
import json
import re
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

API = "https://openebench.bsc.es/api/scientific/staged/Dataset"
EVENT_2022 = {"OEBX00200000" + x for x in
              ["15", "16", "17", "18", "19", "1A", "1B", "1C", "1D", "1E", "1F", "1G", "1H"]}
NAME_RE = re.compile(r"Metric '([^']+)' in challenge '([^']+)' applied to participant '(.+)'")


def get(url):
    err = None
    for _ in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                return fh.read().decode()
        except Exception as e:  # noqa: BLE001
            err = e
    raise RuntimeError(f"failed {url}: {err}")


def main(out):
    ids = [l.strip() for l in get(API).splitlines() if l.startswith("OEBD002")]

    def fetch(i):
        try:
            return json.loads(get(f"{API}/{i}"))
        except Exception:  # noqa: BLE001
            return None

    with ThreadPoolExecutor(16) as ex:
        recs = [r for r in ex.map(fetch, ids) if r]
    n = 0
    with open(out, "w") as fo:
        fo.write("participant\tchallenge\tmetric\tvalue\terror\tdataset_id\n")
        for r in recs:
            if r.get("type") != "assessment" or not set(r.get("challenge_ids", [])) & EVENT_2022:
                continue
            m = NAME_RE.match(r.get("name", ""))
            if not m:
                continue
            met, ch, part = m.groups()
            d = r.get("datalink", {}).get("inline_data", {})
            fo.write(f"{part}\t{ch}\t{met}\t{d.get('value')}\t{d.get('error')}\t{r['_id']}\n")
            n += 1
    print(f"wrote {n} rows from {len(recs)} datasets to {out}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "qfo2022_public_results.tsv")
