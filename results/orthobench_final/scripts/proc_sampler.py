#!/usr/bin/env python
"""Sample cumulative CPU (utime+stime+cutime+cstime of the root process + live descendants) and summed RSS of a
process tree every DT seconds; append TSV lines: epoch, cpu_s, rss_tree_mb, n_procs, load1.
    proc_sampler.py PID OUT.tsv [DT]
Stops when PID exits."""
import os, sys, time
pid, out = int(sys.argv[1]), sys.argv[2]
dt = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
HZ = os.sysconf("SC_CLK_TCK"); PG = os.sysconf("SC_PAGE_SIZE")

def tree(root):
    kids = {}
    for p in os.listdir("/proc"):
        if not p.isdigit():
            continue
        try:
            with open(f"/proc/{p}/stat") as fh:
                s = fh.read()
            ppid = int(s.rsplit(")", 1)[1].split()[1])
            kids.setdefault(ppid, []).append(int(p))
        except Exception:
            pass
    res, st = [], [root]
    while st:
        x = st.pop(); res.append(x); st.extend(kids.get(x, []))
    return res

with open(out, "a") as fo:
    fo.write("epoch\tcpu_s\trss_tree_mb\tn_procs\tload1\n")
    while os.path.exists(f"/proc/{pid}"):
        cpu = 0.0; rss = 0; n = 0
        for p in tree(pid):
            try:
                with open(f"/proc/{p}/stat") as fh:
                    f = fh.read().rsplit(")", 1)[1].split()
                cpu += (int(f[11]) + int(f[12]) + int(f[13]) + int(f[14])) / HZ
                rss += int(f[21]) * PG; n += 1
            except Exception:
                pass
        fo.write("%.1f\t%.1f\t%.0f\t%d\t%s\n" % (time.time(), cpu, rss / 2**20, n, open("/proc/loadavg").read().split()[0]))
        fo.flush()
        time.sleep(dt)
