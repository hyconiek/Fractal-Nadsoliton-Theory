# CONTROLLED-MICRO-TO-HISTOGRAM-300
## One frozen finite-N FIN contract now generates observable histograms directly from the exact microscopic generator

Date: 2026-09-27

Status:
- exact finite-state calculation at N=6, g=5.145228719489142;
- histograms regenerated from the exact leave-one-out Gibbs count generator, not copied from earlier correlation tables;
- exact D12-equivariant basin construction reproduced from the full V_g descent rule;
- finite-N result only; no apparatus assignment and no N->infinity theorem.

## 1. Frozen contract

The contract is:

- N=6 labelled-copy Gibbs system, represented exactly by occupation counts;
- g=5.145228719489142;
- exact leave-one-out heat-bath generator;
- basin map: theta0=(g/N) X7^T n, full V_g descent, then exact D12 orbit propagation;
- readout before compression: 12 localized labels J=0,...,11 plus a residual class R;
- seven-bin Y=cos(2 pi J/12) compression is applied only after checking reflection symmetry.

The state space contains 12,376 count states.

The reconstructed basin partition gives

    localized total mass = 0.992496493817772
    mass per localized basin = 0.082708041151481 (up to roundoff)
    residual stationary mass = 0.007503506182228

This exactly reproduces the accepted N=6 localized-basin mass.

The exact Gibbs stationary law satisfies

    ||pi Q||_inf = 5.28e-17.

So the new observable predictions start from the accepted microscopic law, not from the effective Q12 tables.

## 2. Three microscopic preparations with the same macro label J=0

Three preparations were frozen:

1. DEEP SEED:
   all six copies have label 0, n=(6,0,...,0).

2. EQUILIBRIUM-J0:
   exact Gibbs equilibrium conditioned on the deterministic basin J=0.

3. FLAT-J0:
   uniform distribution over all count states classified as J=0.

All three give exactly the same coarse observation J=0 at t=0.

They are not dynamically equivalent.

| t | TV(seed, equilibrium-J0) | TV(flat-J0, equilibrium-J0) |
|---:|---:|---:|
| 0.5 | 0.03215 | 0.43074 |
| 1 | 0.04517 | 0.50264 |
| 2 | 0.05672 | 0.53203 |
| 4 | 0.06066 | 0.52116 |
| 8 | 0.05568 | 0.47037 |
| 24 | 0.03583 | 0.30448 |
| 64 | 0.01219 | 0.10774 |

Therefore the statement

    "prepare J=0"

is not a complete preparation protocol.

This is a direct PASS of the preparation kill-test proposed after report 299: two microscopic distributions compatible with the same macro label can give measurably different future data.

## 3. Why equilibrium-J0 is the correct reference for the existing projection lane

For equilibrium-J0, the directly generated k=2 response is

    C2(0.5)=0.956168458169
    C2(1)  =0.933060623019
    C2(2)  =0.896603591891
    C2(4)  =0.834172374725
    C2(8)  =0.725021835418

These are the same exact N=6 projected correlations used in report 298.

So report 298 is now embedded in a complete operational chain:

    exact microscopic generator
      -> explicit microscopic preparation
      -> deterministic basin readout
      -> observable histogram / Fourier response.

It is no longer only a post-processing analysis of stored correlations.

## 4. Direct microscopic histogram versus the effective Q12 chain

For the equilibrium-J0 preparation, compare the exact microscopic 13-bin histogram with the memory-renormalized Q12 chain started at J=0.

| t | residual mass | TV full micro vs Q12+R=0 | TV after conditioning on localized states |
|---:|---:|---:|---:|
| 0.5 | 0.006338 | 0.02054 | 0.01437 |
| 1 | 0.007077 | 0.02496 | 0.01819 |
| 2 | 0.007361 | 0.02733 | 0.02048 |
| 4 | 0.007463 | 0.02716 | 0.02056 |
| 8 | 0.007498 | 0.02454 | 0.01852 |
| 16 | 0.007503 | 0.01973 | 0.01475 |
| 24 | 0.007504 | 0.01590 | 0.01176 |
| 32 | 0.007504 | 0.01286 | 0.00939 |
| 64 | 0.007504 | 0.00786 | 0.00389 |

Thus the effective chain is not exact from t=0, but after the memory layer its conditional localized histogram approaches the microscopic prediction to below 0.4% TV on the tested t=64 point.

At the report-297 operating time

    rho = 0.02261891215447
    tau*=rho t*=0.5427059873
    t*=23.9934610291,

we obtain

    TV(micro equilibrium-J0 | localized, Q12) = 0.0117638.

This gives a direct reduction-error budget at the actual fingerprint time.

## 5. Reflection check before seven-bin compression

For the frozen J=0 preparations, the exact full-label histograms satisfy

    p_j(t)=p_{-j}(t)

within numerical roundoff (maximum asymmetry at t* below 7e-16 in the replay).

Therefore the seven-bin Y compression is legitimate for this controlled contract.

This check must be repeated, not assumed, if the preparation or physical readout changes.

## 6. Scientific verdict

Report 300 is positive but changes the interpretation of the effective theory.

The controlled statement is now:

    microscopic law
      + microscopic preparation map
      + basin/readout map
      -> predicted observable histogram.

The macro state J alone does not determine the future process.

For this lane, equilibrium-J0 is the preparation naturally matched to the Mori-Zwanzig / equilibrium-correlation reduction. A deep-seed preparation is operationally simpler but requires its own preparation/slip map. A flat-J0 preparation is not interchangeable with either.

This is a finite-system effective-theory result, not a fundamental source law for A7, g, the update rule or physical space.
