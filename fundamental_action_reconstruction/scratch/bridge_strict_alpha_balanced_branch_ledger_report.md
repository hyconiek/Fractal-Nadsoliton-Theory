# Scratch strict-alpha balanced branch ledger probe

Status: exact balanced-ledger witness for eta=9/5; no selector derivation theorem.

- Weakened rule: replace explicit exponent alphabet `{1,2}` with positive integer exponents plus a balance/equipartition selector over `5` branches summing to `8`.
- Scan size: `35` labelled positive ledgers and `3` canonical ledgers.
- Unique variance and max-gap minimizer: `[2, 2, 2, 1, 1]` with variance `6/25`.
- Exact replay: `(4/3)^3 * (2/3)^2 = 2^8/3^5 = 256/243`, correction residual `0.000e+00`, eta residual `0.000e+00`.
- No false pass: branch count, denominator-3 normalization, and balance selector are still premises; no QW-2191 discharge and no ToE closure.
