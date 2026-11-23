#!/usr/bin/env bash
# submit.sh L-x : verify the step file (rules, not merged upstream, compiles, std axioms), open the PR via the ../submit worktree, tell the desk.
set -euo pipefail
id=$1; R="/Users/Ethan/Library/Application Support/Claude/scratch-workspaces/b53eb492-782f-4b89-a4cd-78a811f342fa/3a48b96d-1945-44a2-a14d-15afce52a49c/scratch-2026-11-23-4f47ca/pequalnp-proofs"
cd "$R"; J=.pnp/$id.json
rel=$(python3 -c "import json;print(json.load(open('$J'))['gap']['file'])"); mod=$(python3 -c "import json;print(json.load(open('$J'))['gap']['module'])")
lease=$(python3 -c "import json;print(json.load(open('$J'))['lease']['id'])")
if grep -qE 'sorry|native_decide|"|^ *#' "$rel"; then echo "RULE VIOLATION"; exit 1; fi
git fetch -q upstream
if git ls-tree upstream/main "$rel" | grep -q .; then echo "ALREADY MERGED UPSTREAM"; exit 1; fi
T=$(mktemp -d); cp "$rel" "$T/c.lean"; echo "#print axioms $mod" >> "$T/c.lean"
ax=$(cd lean && PATH="$HOME/.elan/bin:$PATH" lake env lean "$T/c.lean" 2>&1 | grep -E "error|depends on axioms" || true)
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
