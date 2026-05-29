# Scratch strict-alpha binary/ternary ledger probe

Status: exact branch-ledger witness for eta=9/5; no branch-law derivation theorem.

- Exact replay: `q^5=1.053497942386831=256/243=2^8/3^5`.
- Candidate finite branch law: five branches, denominator `3` per branch, and binary exponents only `1` or `2`.
- Unique ledger under that rule: `[2, 2, 2, 1, 1]`, i.e. `(4/3)^3*(2/3)^2`, with correction residual `0.000e+00` and eta residual `0.000e+00`.
- Looser scan warning: `11` canonical ledgers exist for length five, sum eight, entries 0..4; uniqueness requires the extra one-or-two-bit branch rule.
- No false pass: no kernel identity, no legacy role transfer, no strict branch-law theorem, no QW-2191 discharge, no ToE closure.
