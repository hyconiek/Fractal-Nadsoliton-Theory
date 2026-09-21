# MP7-040 — Fixed-fixture phase roots versus full equilibria

Status: **PROVED_INTERVAL_ASSISTED within the current MP7 working campaign; serialized at handoff time.**

Scope: the exact-decimal fixed-amplitude fixture used by the accepted 60-root full-phase census. This statement does **not** vary the amplitudes and does not classify nearby amplitude-phase stationary points.

## Claim

None of the exactly 60 certified phase critical points at the fixed fixture is a full equilibrium satisfying

`theta = g * mu(theta)`

for any positive scalar `g`.

Equivalently, the previously certified 60 phase roots are phase-critical points at frozen amplitudes; they are not 60 full stationary states of the rank-seven mediator model.

## Split of the proof

1. **54 generic roots.** On each root's certified full-phase uniqueness collar, an interval enclosure of the projected stationarity residual

   `P_{theta^perp} mu`

   excludes zero. The weakest separating component still has an absolute lower separation of approximately `8.89e-4`.

2. **6 symmetry-locked roots.** These are the odd-translation orbit of an aligned field with the alternating amplitude sign flipped. Their root IDs in the R7N full-phase catalog are

   `11, 14, 18, 40, 45, 49`.

   Representative phase locks are the odd label translations of the aligned configuration; for example root 11 is at `(pi/2,0,3pi/2)` up to certified root-box error. After translating to the aligned positive representative, full stationarity would require the four radial ratios `rho_i = mu_i/theta_i = 1/g` to coincide. Outward interval evaluation separates them. The tightest comparison was between modes 5 and 6:

   `rho_5 - rho_6 in [1.00428e-10, 1.00446e-10] > 0`.

   Hence no common positive `g` exists.

## Inputs

- `R7N-041_full_uniqueness_collars.json`: 60 full-phase roots with certified collars.
- Exact-decimal fixed amplitudes inherited from the R7N phase fixture.
- The audited rank-seven feature/softmax model used throughout R7N/R7O3/MP7.

## Nonconclusions

- This does not say that the 60 phase roots disappear when amplitudes are allowed to move.
- It does not exclude nearby full amplitude-phase equilibria.
- It does not classify all stationary points of the full seven-dimensional model.
- It does not change the exact 60-point phase census at the frozen fixture.
- A local amplitude-to-phase continuation (MP7-041) is the correct next tool if one wants to follow one of these phase branches as amplitudes vary.
