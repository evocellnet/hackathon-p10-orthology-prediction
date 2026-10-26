#!/usr/bin/env python3
"""Download the public QfO benchmark results (all participants, all challenges) from the
OpenEBench scientific GraphQL API and write them as raw JSON + a tidy TSV.

Events: QfO 2022 = OEBE0020000003, 2020 = OEBE0020000002, 2018 = OEBE0020000001, 2011 = OEBE0020000000
Usage: fetch_openebench_results.py [--event OEBE0020000003] [--outdir reference_results]
"""
import argparse, csv, datetime, json, os, urllib.request

API = "https://openebench.bsc.es/api/scientific/graphql"

def gql(q):
    req = urllib.request.Request(API, data=json.dumps({"query": q}).encode(),
                                 headers={"Content-Type": "application/json"})
    return json.load(urllib.request.urlopen(req, timeout=120))["data"]

p = argparse.ArgumentParser()
p.add_argument("--event", default="OEBE0020000003")
p.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__), "..", "reference_results"))
a = p.parse_args()
os.makedirs(a.outdir, exist_ok=True)
ev = gql('{ getBenchmarkingEvents(benchmarkingEventFilters:{id:"%s"}) { _id name orig_id } }' % a.event)
evname = ev["getBenchmarkingEvents"][0]["orig_id"].replace(":", "_")
chals = gql('{ getChallenges(challengeFilters:{benchmarking_event_id:"%s"}) { _id acronym orig_id } }' % a.event)["getChallenges"]
raw, rows = {"event": ev, "retrieved": datetime.datetime.now().isoformat(), "api": API, "challenges": {}}, []
for c in sorted(chals, key=lambda c: c["acronym"]):
    ds = gql('{ getDatasets(datasetFilters:{challenge_id:"%s", type:"aggregation"}) { _id orig_id datalink { inline_data } } }' % c["_id"])["getDatasets"]
    raw["challenges"][c["acronym"]] = {"challenge": c, "datasets": []}
    for d in ds:
        inl = d["datalink"]["inline_data"]
        inl = json.loads(inl) if isinstance(inl, str) else inl
        raw["challenges"][c["acronym"]]["datasets"].append({"_id": d["_id"], "orig_id": d["orig_id"], "inline_data": inl})
        vis = inl.get("visualization", {})
        for part in inl.get("challenge_participants", []):
            rows.append({"event": evname, "challenge": c["acronym"], "dataset": d["orig_id"],
                         "x_metric": vis.get("x_axis", ""), "y_metric": vis.get("y_axis", ""),
                         "tool": part.get("tool_id") or part.get("participant_id"),
                         "x": part.get("metric_x"), "stderr_x": part.get("stderr_x"),
                         "y": part.get("metric_y"), "stderr_y": part.get("stderr_y"),
                         "optimization": vis.get("optimization", "")})
    print(c["acronym"], len(ds), "datasets")
base = os.path.join(a.outdir, f"openebench_{evname}")
json.dump(raw, open(base + "_raw.json", "w"), indent=1)
with open(base + ".tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
print("wrote", base + ".tsv", len(rows), "rows")
