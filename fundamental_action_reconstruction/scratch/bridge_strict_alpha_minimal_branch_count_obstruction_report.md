# Scratch strict-alpha minimal branch-count obstruction probe

Status: irreducible branch-count certificate for eta=9/5; no physical branch theorem.

- Model rule: `q_model=(2^n/3^m)^(1/m)=2^(n/m)/3`; exact target requires `n/m=8/5`, equivalently `5n=8m`.
- General solution: `m=5*t`, `n=8*t`; bounded scan found `8` exact matches through `m<=40`, `n<=80`.
- Minimal exact match: `m=5`, `n=8`, eta residual `0.000e+00`.
- Honest read: five branches are arithmetically irreducible inside this branch model, but the model itself and the balanced-ledger selector are still not derived.
- No false pass: no denominator-3 branch theorem, no selector discharge, no QW-2191 discharge, no ToE closure.
