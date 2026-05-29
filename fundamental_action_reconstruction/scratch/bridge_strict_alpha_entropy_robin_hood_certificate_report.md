# Scratch strict-alpha entropy Robin-Hood certificate probe

Status: constructive fixed-labelled entropy certificate for eta=9/5; no strict selector theorem.

- Weight law: `W(e)=8!/prod_i(e_i!)`; a balancing move `a -> a-1`, `b -> b+1` has exact ratio `a/(b+1)`.
- Certificate path: `(4,1,1,1,1) -> (3,2,1,1,1) -> (2,2,2,1,1)`, with labelled weights `1680 -> 3360 -> 5040`.
- Terminal ledger: `[2, 2, 2, 1, 1]`, entropy `1.559581156260`, eta residual `0.000e+00`.
- Honest read: fixed-labelled multiplicity constructively forces the balanced ledger, but fixed-labelled strict channels are still a missing selector premise.
- No false pass: no fixed-labelled strict selector theorem, no QW-2191 discharge, no ToE closure.
