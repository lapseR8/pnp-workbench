#!/usr/bin/env bash
# release.sh L-x reason : release a step and delete its skeleton file.
R="/Users/Ethan/Library/Application Support/Claude/scratch-workspaces/b53eb492-782f-4b89-a4cd-78a811f342fa/3a48b96d-1945-44a2-a14d-15afce52a49c/scratch-2026-11-23-4f47ca/pequalnp-proofs"
cd "$R" && tools/pnp.sh release "$1" "$2" && rm -f "$(python3 -c "import json;print(json.load(open('.pnp/$1.json'))['gap']['file'])")" && rm -f ".pnp/$1.json"
