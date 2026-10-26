#!/usr/bin/env bash
# Run one competitor command under /usr/bin/time -v with node-load logging (uptime every 60 s).
# Usage: run_competitor_timed.sh <label> <outdir-for-logs> -- <command...>
# Writes <outdir>/time.txt (time -v), <outdir>/load.log, <outdir>/cmd.log, <outdir>/start_end.txt
set -u
label="$1"; od="$2"; shift 3
mkdir -p "$od"
( while true; do echo "$(date -Iseconds) $(uptime)"; sleep 60; done ) > "$od/load.log" 2>&1 &
LP=$!
date -Iseconds > "$od/start.txt"
/usr/bin/time -v -o "$od/time.txt" "$@" > "$od/cmd.log" 2>&1
rc=$?
date -Iseconds > "$od/end.txt"
echo "$rc" > "$od/exit_code.txt"
kill $LP 2>/dev/null
exit $rc
