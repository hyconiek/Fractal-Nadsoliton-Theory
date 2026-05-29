# Scratch strict-alpha selector orbit-weight threshold probe

Status: exact selector-threshold discriminator for eta=9/5; no strict selector theorem.

- Selector family: `Score_gamma(e)=log W(e)+gamma*log O(e)`, interpolating fixed-labelled (`gamma=0`) and orbit-aggregate (`gamma=1`) conventions.
- Critical threshold: `gamma_c=log(3/2)/log(2)=0.584962500721`.
- Winners: below threshold `[2, 2, 2, 1, 1]`; above threshold `[3, 2, 1, 1, 1]`; at `gamma=1`, `[3, 2, 1, 1, 1]`.
- Honest read: the balanced eta-ledger is stable only for a below-threshold selector convention; this is a quantitative selector proof obligation, not closure.
- No false pass: no strict derivation of gamma=0 or gamma<gamma_c, no QW-2191 discharge, no ToE closure.
