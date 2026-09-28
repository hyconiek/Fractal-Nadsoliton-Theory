# JOINT-FINGERPRINT-PROTOCOL-302
## A single frozen N=6 protocol separates preparation, short memory, late rate geometry and the same-rho comparator

Date: 2026-09-27

Status:
- exact finite-N microscopic FIN predictions at N=6;
- same-rho/same-total-exit comparator from report 293/297;
- explicit 5% symmetric categorical readout model;
- local timing/readout/g sensitivity checks;
- not a laboratory validation and not a universal minimax test against all alternatives.

## 1. Calibration/test split

The protocol is deliberately split.

CALIBRATION:
- freeze N=6, g=5.145228719489142;
- use the exact equilibrium-J0 microscopic preparation;
- retain all 12 localized labels plus the residual class;
- verify reflection symmetry before compression;
- calibrate the slow clock rho from an independent late-time data set;
- calibrate the readout confusion model independently.

TEST 1 — EARLY MEMORY/PREPARATION:
- use multiple early times;
- test variation of the k=2 logarithmic slope.

TEST 2 — LATE FIN FINGERPRINT:
- use tau=rho t=0.5427059873;
- compress to the seven-bin Y=cos(2 pi J/12) histogram only after the symmetry check;
- compare with the predeclared same-rho/same-total-exit comparator.

No parameter is re-fitted separately at each test time.

## 2. Early test from the exact microscopic histogram

For the equilibrium-J0 preparation:

    C2(0.5)=0.956168458169
    C2(1)  =0.933060623019
    C2(2)  =0.896603591891.

The logarithmic slopes are

    s(0.5,1) = -0.04892786765
    s(1,2)   = -0.03985633734.

Thus

    Delta_s = s(1,2)-s(0.5,1)
            = 0.00907153031.

A pure single exponential has Delta_s=0.

Normalize by the independently calibrated slow clock:

    Delta_s/rho = 0.4011.

This dimensionless ratio is invariant under one common multiplicative rescaling of the time unit.

For the symmetric time-independent categorical readout model used in report 297, every nonzero Fourier mode is multiplied by one constant factor. That constant cancels from logarithmic slope differences. Thus the mean early-memory statistic is exactly insensitive to this declared readout attenuation; sampling variance still increases.

## 3. Late seven-bin fingerprint generated from the microscopic law

The operating point is

    rho = 0.02261891215447
    tau* = 0.5427059873
    t* = 23.9934610291.

With 5% symmetric readout error, the direct microscopic equilibrium-J0 FIN prediction has TV distance

    0.0740883

from the report-293 same-rho/same-total-exit comparator.

Chernoff information is

    C = 0.007307535614.

For equal priors, the standard bound

    Pe <= 0.5 exp(-M C)

falls below 5% at

    M = 316 independent shots.

This number is specific to this declared comparator, preparation and noise model.

The memory-renormalized effective Q12 prediction gives C=0.00727880, very close to the direct microscopic value, but the microscopic histogram is now the primary prediction and Q12 is only its controlled reduction.

## 4. Explicit reduction-error budget at the test time

At t* the exact microscopic equilibrium-J0 histogram, after conditioning on localized outcomes, differs from Q12 by

    TV = 0.0117638.

The direct FIN-versus-comparator seven-bin separation is

    TV = 0.0740883.

Thus at this point the reduction discrepancy is about 15.9% of the model-separation TV scale.

It is not negligible enough to omit from a precision budget, but it is substantially smaller than the declared comparator separation.

## 5. Readout robustness

For the same frozen microscopic prediction and comparator:

| symmetric error eta | Chernoff C | 5% Chernoff-bound shots |
|---:|---:|---:|
| 0 | 0.0087331 | 264 |
| 0.02 | 0.0081314 | 284 |
| 0.05 | 0.0073075 | 316 |
| 0.10 | 0.0061124 | 377 |
| 0.20 | 0.0042366 | 544 |

So this specific fingerprint survives substantial symmetric readout dilution.
Asymmetric or time-dependent readout is NOT covered by this cancellation and must be calibrated separately.

## 6. Timing robustness

Keeping the protocol otherwise fixed and shifting the late test time by +/-20% gives:

| t/t* | Chernoff C | 5% bound shots |
|---:|---:|---:|
| 0.8 | 0.0072046 | 320 |
| 0.9 | 0.0073025 | 316 |
| 1.0 | 0.0073075 | 316 |
| 1.1 | 0.0072374 | 319 |
| 1.2 | 0.0071072 | 324 |

The late protocol is therefore broad rather than sharply tuned.

## 7. Local microscopic-parameter sensitivity

As a local sensitivity check only, keep the nominal basin readout fixed and recompute the exact microscopic generator at

    g0-0.01 and g0+0.01.

At t* the seven-bin prediction moves by about

    TV(nominal, g0-0.01) = 0.00511
    TV(nominal, g0+0.01) = 0.00508.

The symmetric finite-difference local sensitivity is about

    0.510 TV units per unit g.

This is NOT a validated experimental uncertainty range and must not be extrapolated far from g0. It merely provides the first derivative-scale entry in the error ledger.

## 8. Preparation uncertainty is currently the dominant unresolved nuisance

At t*, with the same 5% readout model:

    TV(seed, equilibrium-J0) = 0.03414
    TV(flat-J0, equilibrium-J0) = 0.29010
    TV(equilibrium-J0 FIN, comparator) = 0.07409.

Thus even the deep-seed versus equilibrium preparation ambiguity is about 46% of the declared model-separation TV scale.
The flat-J0 ambiguity is almost four times larger than the model signal.

Therefore a protocol that says only "prepare J=0" is not falsifiable at the claimed precision.
The microscopic preparation map must be part of the model contract.

## 9. Verdict

The combined 297-298 idea survives a stricter direct-microscopic implementation PROVIDED the preparation is frozen.

The two tests probe different things:

    early multi-time response -> unresolved memory / preparation layer;
    late seven-bin histogram -> FIN-specific shell-rate geometry beyond rho.

The protocol is robust to the declared symmetric readout noise and broad timing shifts, and it now has an explicit finite-N reduction-error entry.

The dominant remaining uncertainty is not the late clock. It is the microscopic preparation map and, after that, the unsourced empirical assignment of N, g and labels to a real system.

## 10. Next research

The highest-value next step is a true held-out prediction:

    303 — HELD-OUT-N-MICROSCOPIC-PREDICTION

Freeze the procedure here, move to N=7 (then N=8), regenerate the microscopic histogram directly, and compare it with the already-existing q_d(N) effective prediction WITHOUT fitting the test histogram.

A second parallel target is:

    PREPARATION-MAP-CERTIFICATE

derive or bound the map from an operationally simple deep-seed preparation to the equilibrium/slip-corrected effective initial condition.
