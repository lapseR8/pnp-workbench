#!/usr/bin/env bash
SELF="$(cd "$(dirname "$0")" && pwd)"
# submit.sh L-x : verify the step file (rules, not merged upstream, compiles, std axioms), open the PR via the ../submit worktree, tell the desk.
set -euo pipefail
id=$1; R="/Users/Ethan/Library/Application Support/Claude/scratch-workspaces/b53eb492-782f-4b89-a4cd-78a811f342fa/3a48b96d-1945-44a2-a14d-15afce52a49c/scratch-2026-11-23-4f47ca/pequalnp-proofs"
cd "$R"; J=.pnp/$id.json
rel=$(python3 -c "import json;print(json.load(open('$J'))['gap']['file'])"); mod=$(python3 -c "import json;print(json.load(open('$J'))['gap']['module'])")
lease=$(python3 -c "import json;print(json.load(open('$J'))['lease']['id'])")
if grep -qE 'sorry|native_decide|"|^ *#' "$rel"; then echo "RULE VIOLATION"; exit 1; fi
git fetch -q upstream 2>/dev/null || echo "note: git fetch upstream failed; using API"
if gh api "repos/brilliant-hq/pequalnp-proofs/contents/$rel" --silent 2>/dev/null; then echo "ALREADY MERGED UPSTREAM"; exit 1; fi
T0=$(mktemp -d); curl -s https://data.pequalnp.com/data/desk/gaps.json -o "$T0/gaps.json"
python3 - "$T0/gaps.json" "$id" "$T0/st.txt" <<'PY' || exit 1
import json, sys
d = json.load(open(sys.argv[1])); i = sys.argv[2]
gaps = d["gaps"]; gaps = gaps if isinstance(gaps, list) else list(gaps.values())
def ids(x): return {(g.get("id") if isinstance(g, dict) else g) for g in (x if isinstance(x, list) else list(x))}
g = next((g for g in gaps if g.get("id") == i), None)
if g is None or i in ids(d.get("withdrawn", [])) or i in ids(d.get("checked", [])):
    print("GUARD: step not posted (withdrawn/checked/missing)"); sys.exit(1)
if (g.get("weight") or 0) <= 1:
    print("GUARD: weight 1 (stood down)"); sys.exit(1)
open(sys.argv[3], "w").write(g["statement"])
print("GUARD: posted, weight", g.get("weight"))
PY
c=$(curl -s -o /dev/null -w '%{http_code}' -X POST "https://pequalnp.com/desk/beat?lease=$lease")
[ "$c" = 200 ] || { echo "GUARD: lease $c (not held)"; exit 1; }
gh pr list -R brilliant-hq/pequalnp-proofs --search "$id in:title" --state merged --json number --jq '.[].number' | grep -q . && { echo "GUARD: a PR for $id is already merged"; exit 1; }
echo "running official checker..."
node tools/check-lean.mjs "$rel" "$T0/st.txt" check > "$T0/check.txt" 2>&1; rc=$?
tail -2 "$T0/check.txt"; [ $rc = 0 ] || { echo "GUARD: official checker rejected (exit $rc)"; exit 1; }
T=$(mktemp -d); cp "$rel" "$T/c.lean"; echo "#print axioms $mod" >> "$T/c.lean"
if [ -S /tmp/claude-502/leanserver.sock ]; then
  ax=$(python3 "$SELF/leanq.py" "$rel" --axioms "$mod" --all 2>&1 | grep -E "^error|depends on axioms|SERVER ERROR" || true)
else
  ax=$(cd lean && lake env lean "$T/c.lean" 2>&1 | grep -E "error|depends on axioms" || true)
fi
echo "$ax"; echo "$ax" | grep -q error && { echo "COMPILE ERROR"; exit 1; }
echo "$ax" | grep -q "depends on axioms: \[propext, Classical.choice, Quot.sound\]" || { echo "AXIOM PROBLEM"; exit 1; }
W="$R/../submit"; [ -d "$W" ] || git worktree add -q "$W" upstream/main
git -C "$W" checkout -q --detach upstream/main; git -C "$W" checkout -q -B "$id"
mkdir -p "$(dirname "$W/$rel")"; cp "$rel" "$W/$rel"; git -C "$W" add "$rel"
git -C "$W" -c user.name=lapseR8 -c user.email=269838972+lapseR8@users.noreply.github.com commit -q -m "$id" -m "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
git -C "$W" -c credential.helper= -c credential.helper='!gh auth git-credential' push -q -f -u origin "$id" 2>&1 | grep -v '^remote:' || true
gh pr create --repo brilliant-hq/pequalnp-proofs --head "lapseR8:$id" --title "$id" --body "**X handle**: @lapseR8wx
**LinkedIn profile**: https://www.linkedin.com/in/ethan-roberts-production/

🤖 Generated with [Claude Code](https://claude.com/claude-code)"
curl -s -X POST "https://pequalnp.com/desk/release?lease=$lease&reason=submitted"; echo
