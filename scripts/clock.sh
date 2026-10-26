#!/usr/bin/env bash
# Print elapsed and remaining time vs START.txt (48 h budget).
root="$(cd "$(dirname "$0")/.." && pwd)"
start=$(head -1 "$root/START.txt"); s=$(date -d "$start" +%s); n=$(date +%s)
e=$(( n - s )); r=$(( s + 172800 - n ))
printf "now %s | elapsed T+%dh%02dm | remaining %dh%02dm\n" "$(date -Iseconds)" $((e/3600)) $((e%3600/60)) $((r/3600)) $((r%3600/60))
