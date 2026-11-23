# pnp-workbench: playbook for continuing the pequalnp.com effort

Owner GitHub handle: lapseR8. Goal: as many ACCEPTED proofs in brilliant-hq/pequalnp-proofs as possible.
Rules: https://pequalnp.com/join.txt (read it). Fork: github.com/lapseR8/pequalnp-proofs.

## Approved by the owner (standing)
- Claim steps at the desk, open PRs to brilliant-hq/pequalnp-proofs under lapseR8.
- Commit as `lapseR8 <269838972+lapseR8@users.noreply.github.com>` (never the personal email).
- PR description, exactly (only these lines are read for credit):
  **X handle**: @lapseR8wx
  **LinkedIn profile**: https://www.linkedin.com/in/ethan-roberts-production/
- Heavy compute goes to GitHub Actions in this repo (public, so free minutes), not the owner's laptop.
- Token usage is not the constraint; accepted output is.

## Loop
1. Setup (cloud or fresh machine): `gh repo clone lapseR8/pequalnp-proofs && cd pequalnp-proofs && git remote add upstream https://github.com/brilliant-hq/pequalnp-proofs || true; git fetch upstream`; install elan; `cd lean && lake exe cache get` (~10 GB); `lake build PEqualNP.Defs` (and any gap.imports).
2. Claim: `tools/pnp.sh next --another` (up to 5 held). Read STEP / KNOWN / SKETCH / RANGE.
3. Triage in under a minute:
   - Short with a sketch or reusable imports (like L-meta-748, L-meta-733): PROVE IT YOURSELF IMMEDIATELY. Don't delegate; others race the same step and win by minutes.
   - Known route but long (TRACTABLE rows with a paper route): chunk it. A planner writes a compiling skeleton of `sorry` lemmas, parallel workers fill lemmas, an assembler merges. Others won 536/540/215 with 7-10k line files after we quit in minutes. Persist; don't quit early.
   - Finite search: SAT via `.github/workflows/batch.yml` (input spec `TAG K N GRAPHJSON;...`, K=99 means no OR limit), then verify the witness in Lean with `decide` (8 inputs ~50 s).
   - Open conjecture with no finite handle: release with an honest reason.
4. Lease: renew every 10 min (`tools/pnp.sh beat <step>`), or put lease IDs in secret `LEASES` and run workflow `leases`. A 410 means stop at once (settled elsewhere).
5. Verify: `cd lean && lake env lean PEqualNP/Outside/L_x.lean`; check axioms on a temp copy with `#print axioms` (only propext, Classical.choice, Quot.sound); no sorry/native_decide/double quotes/# commands; helpers inside `namespace PEqualNP.Outside.L_x`; theorem last. Check the step is not already in `upstream/main`.
6. Submit: `git checkout -b L-x upstream/main`, add the one file, commit `L-x`, push to origin, `gh pr create --repo brilliant-hq/pequalnp-proofs --head lapseR8:L-x --title L-x --body <credit lines>`, then `curl -X POST 'https://pequalnp.com/desk/release?lease=<id>&reason=submitted'`.
7. Log every outcome in README.md (results table) and save proofs under `proofs/`.

## State at handoff (2026-11-23 ~17:45 UTC)
- PR brilliant-hq/pequalnp-proofs#1170 (L-meta-748) open, awaiting checker verdict (comment on the PR).
- Held: L-gct-179 (find a 7-dim complex space of 4x5 matrices where sum over injections of prod(1+t) == 120). Lease renewed by workflow `leases` until ~22:00 UTC. A numeric search agent was running locally and stops with the laptop. Hint: study upstream L_gct_215 (proves the bound 7) for which case could be tight.
- Lost races: 733 (our proof in proofs/), 719, 586, 583, 619. Others proved 536/540/215, which we had released.

## Update 2026-11-23 19:40 UTC
- All our PRs so far lost races: #1170 (748), #1831 (948), #1969 (239). Winners PR within minutes of a row going live (fededevi, stierms, bz, ...). The desk even hands out steps already merged.
- RULE: right after claiming, run `gh pr list -R brilliant-hq/pequalnp-proofs --search "<STEP> in:title" --state all`. Invest only if there is no earlier PR (or all closed), and renew EVERY held lease from GitHub (secret LEASES + workflow `leases`) immediately.
- Explicit-circuit steps ("decide on all 256 inputs"): template `and_intros; all_goals (intro C hC; subst hC; refine ⟨?_, rfl, ?_⟩; simp [FanInTwo]; decide)` with `set_option maxRecDepth 100000 in set_option maxHeartbeats 0 in`, ~80 s for 3 circuits.
- Held: L-gct-179; numeric search workflow `gct179` (agent-built) running.
