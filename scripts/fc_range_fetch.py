#!/usr/bin/env python
"""Build a mini Foldcomp DB (data/index/lookup/dbtype) containing only the requested keys of a remote
Foldcomp DB by HTTP range requests (S3). Usage:
  fc_range_fetch.py <remote_base_url> <full.index> <full.lookup_or_key2name.tsv> <keys.txt> <out_prefix> [--threads 32]
Only needs the (local) .index of the full DB; reads ~5 KB per entry instead of the 1.2 TB data file."""
import sys, urllib.request, concurrent.futures as cf, time
url, idxf, lookf, keysf, out = sys.argv[1:6]
thr = int(sys.argv[sys.argv.index('--threads') + 1]) if '--threads' in sys.argv else 32
keys = [l.strip() for l in open(keysf) if l.strip()]
kset = set(keys)
idx = {}
with open(idxf) as f:
    for l in f:
        k = l.split('\t', 1)[0]
        if k in kset:
            _, o, n = l.split()
            idx[k] = (int(o), int(n))
names = {}
with open(lookf) as f:
    for l in f:
        t = l.rstrip('\n').split('\t')
        if t[0] in kset:
            names[t[0]] = t[1]
def get(k):
    o, n = idx[k]
    for att in range(6):
        try:
            r = urllib.request.Request(url, headers={'Range': f'bytes={o}-{o+n-1}'})
            d = urllib.request.urlopen(r, timeout=60).read()
            if len(d) == n:
                return k, d
        except Exception:
            time.sleep(1 + att)
    return k, None
t0 = time.time(); off = 0; bad = 0
with open(out, 'wb') as df, open(out + '.index', 'w') as ix, open(out + '.lookup', 'w') as lk, cf.ThreadPoolExecutor(thr) as ex:
    for i, (k, d) in enumerate(ex.map(get, keys)):
        if d is None:
            bad += 1; continue
        df.write(d); ix.write(f'{k}\t{off}\t{len(d)}\n'); lk.write(f'{k}\t{names[k]}\t0\n'); off += len(d)
        if i % 20000 == 0: print(i, f'{time.time()-t0:.0f}s', file=sys.stderr)
open(out + '.dbtype', 'wb').write(b'\x0c\x00\x00\x00')
print('fetched', len(keys) - bad, 'failed', bad, f'{time.time()-t0:.0f}s', file=sys.stderr)
