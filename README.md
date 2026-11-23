# pnp-workbench

Working tools and results for contributions to brilliant-hq/pequalnp-proofs (pequalnp.com).

## Tools
- `sat/gen.py N K GRAPH out.cnf out.map` writes a CNF: an exact N-gate fan-in-2 AND/OR circuit
  with at most K OR gates computing the edge function of GRAPH (JSON `{"n":..,"E":[[u,v],..]}`).
  Strict operands, every gate used, no duplicate gates. Solve with `kissat -q out.cnf > out.out`.
- `sat/decode.py out.map out.out` prints the circuit and checks it on every input.
- `sat/search2.py K lo hi GRAPH` is the same search via python-sat/CaDiCaL (allows copy gates, so
  UNSAT at N also rules out every smaller size).
- `lean-tests/T8.lean` shows that `decide` checks an 8-input circuit on all 256 inputs (~50 s).

Setup: `python3 -m venv venv && venv/bin/pip install python-sat`, `brew install kissat`.

## Results log (2026-11-23)
| Step | Result |
|---|---|
| L-meta-619 (wheel W5, <=8 gates) | UNSAT at 8 gates (search2.py, copies allowed, so covers <=8). Cross-check at 9 running. |
| L-meta-586 (two 5-cycles sharing an edge, <=7 OR) | 8 OR found (12 gates). 7 OR: kissat runs at N=11..17 in progress. |
| L-meta-536 | released: missing_lemma (needs uncommitted step lemmas 538-540) |
| L-meta-606 | released: missing_lemma (Mehlhorn-Galil/Pippenger-style lower bound) |
| L-meta-540 | released: out_of_time (two-gate elimination engine needed) |
| L-gct-215 | released: out_of_time (repo methods give <=8 or <=9, not 7) |

8-OR circuit for L-meta-586 (wires 0-7 inputs):
`or(4,7) or(5,7) and(6,9) or(0,2) or(2,4) or(5,11) and(0,8) or(10,14) and(3,12) or(15,16) and(1,13) or(17,18)`
