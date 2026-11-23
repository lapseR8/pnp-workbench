#!/usr/bin/env bash
# claim.sh [--hard] : claim one step, triage it in seconds, and put every live lease on GitHub renewal.
# Run from anywhere. Prints STEP info, upstream/PR status, and a verdict line.
set -uo pipefail
R="/Users/Ethan/Library/Application Support/Claude/scratch-workspaces/b53eb492-782f-4b89-a4cd-78a811f342fa/3a48b96d-1945-44a2-a14d-15afce52a49c/scratch-2026-11-23-4f47ca/pequalnp-proofs"
WB=lapseR8/pnp-workbench
cd "$R" || exit 1
before=$(ls .pnp/*.json 2>/dev/null | sort)
tools/pnp.sh next --another "$@" | grep -E "^(STEP|FILE|KNOWN|SKETCH|RANGE|[0-9]+:|YOUR|LEASE)" | cut -c1-500
after=$(ls .pnp/*.json 2>/dev/null | sort)
new=$(comm -13 <(echo "$before") <(echo "$after") | head -1)
[ -z "$new" ] && new=$(ls -t .pnp/L-*.json | head -1)
id=$(basename "$new" .json)
rel=$(python3 -c "import json;print(json.load(open('$new'))['gap']['file'])")
git fetch -q upstream
if git ls-tree upstream/main "$rel" | grep -q .; then echo "VERDICT: ALREADY MERGED UPSTREAM -> release"; fi
prs=$(gh pr list -R brilliant-hq/pequalnp-proofs --search "$id in:title" --state all --json number,author,createdAt,state --jq 'sort_by(.createdAt)|.[]|"\(.number) \(.author.login) \(.createdAt) \(.state)"')
[ -n "$prs" ] && { echo "COMPETING PRs:"; echo "$prs"; } || echo "NO COMPETING PRs"
# collect all live leases and refresh GitHub renewal
live=""
for f in .pnp/L-*.json; do
  l=$(python3 -c "import json;print(json.load(open('$f'))['lease']['id'])")
  c=$(curl -s -o /dev/null -w '%{http_code}' -X POST "https://pequalnp.com/desk/beat?lease=$l"); sleep 1
  [ "$c" = 200 ] && live="$live $l"
done
printf %s "$live" | gh secret set LEASES -R $WB
for r in $(gh run list -R $WB -w leases --status in_progress --json databaseId --jq '.[].databaseId'); do gh run cancel "$r" -R $WB >/dev/null; done
gh workflow run leases -R $WB >/dev/null && echo "LEASES ON GITHUB:$(echo $live | wc -w | tr -d ' ') live"
