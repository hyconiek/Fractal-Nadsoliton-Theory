# Scratch strict-alpha entropy selector discriminator probe

Status: selector discriminator for eta=9/5; no strict selector theorem.

- Candidate space: positive five-branch ledgers summing to `8`: `(4,1,1,1,1)`, `(3,2,1,1,1)`, `(2,2,2,1,1)`.
- Fixed-labelled entropy winner: `[2, 2, 2, 1, 1]` with labelled multinomial count `5040`.
- Unlabelled orbit-aggregate winner: `[3, 2, 1, 1, 1]` with aggregate count `67200`.
- Honest read: labelled Shannon/multinomial selection supports the balanced ledger, but orbit aggregation selects a different ledger; the selector convention is a real proof obligation.
- No false pass: no fixed-labelled strict selector theorem, no QW-2191 discharge, no ToE closure.
