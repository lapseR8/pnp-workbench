#!/usr/bin/env bash
# run_one.sh TAG K N GRAPH : one SAT search; writes results/TAG_N<N>.txt and pushes it.
set -u
tag=$1 k=$2 n=$3 g=$4
d=$(mktemp -d)
python -I sat/gen.py "$n" "$k" "$g" "$d/k.cnf" "$d/k.map" > /dev/null
start=$(date +%s)
timeout 19800 kissat -q "$d/k.cnf" > "$d/k.out"; rc=$?
secs=$(( $(date +%s) - start ))
mkdir -p results
verdict=$(grep -m1 '^s ' "$d/k.out" || echo "s TIMEOUT (exit $rc)")
f=results/${tag}_N${n}.txt
echo "N=$n K=$k: $verdict, ${secs}s, kissat $(kissat --version)" > "$f"
grep -q '^s SATISFIABLE' "$d/k.out" && python -I sat/decode.py "$d/k.map" "$d/k.out" >> "$f"
cat "$f"
for i in 1 2 3 4 5 6 7 8; do
  ( flock 9
    git add "$f" && git commit -q -m "$tag search result N=$n" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
    git pull -q --rebase && git push -q ) 9>/tmp/git.lock && break
  sleep $((RANDOM % 20 + 5))
done
