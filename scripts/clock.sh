#!/usr/bin/env bash
# Elapsed/remaining budget. Supports pause/resume lines in START.txt ("resume: <iso> budget_remaining=..." + "new deadline <iso>").
root="$(cd "$(dirname "$0")/.." && pwd)"
start=$(head -1 "$root/START.txt"); s=$(date -d "$start" +%s); n=$(date +%s)
dl=$(grep -o 'new deadline [^)]*' "$root/START.txt" | tail -1 | cut -d' ' -f3)
if [ -n "$dl" ]; then d=$(date -d "$dl" +%s); else d=$(( s + 172800 )); fi
r=$(( d - n )); used=$(( 172800 - r ))
printf "now %s | budget used %dh%02dm (pause excluded) | remaining %dh%02dm | deadline %s\n" "$(date -Iseconds)" $((used/3600)) $((used%3600/60)) $((r/3600)) $((r%3600/60)) "$(date -Iseconds -d @$d)"
