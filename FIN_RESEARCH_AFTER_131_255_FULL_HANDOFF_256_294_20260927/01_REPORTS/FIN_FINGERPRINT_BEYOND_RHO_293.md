# FIN-FINGERPRINT-BEYOND-RHO-293
## The resolved 12-state effective dynamics contains dimensionless FIN-specific response data not fixed by the coarse Z3 clock

Date: 2026-09-27

Status:
exact algebra inside the reconstructed 12-state D12-circulant effective chain;
numerical q_d inputs for N=3..8 from report 214.

The exact Z3 quotient fixes

    k_Z3
      =
    q1+q2+q4+q5,

and

    rho
      =
    3 k_Z3.

If an experiment or reduced description observes only the Z3 process, all shell decompositions with the same k_Z3 are indistinguishable at that level.

But the resolved 12-state chain has additional Fourier modes.

## 1. Hidden/resolved-mode eigenvalues

For k=0,...,6,

    lambda_k
      =
    sum_(d=1)^5
    2 q_d[
      cos(2 pi k d/12)-1
    ]
    +
    q6[
      cos(pi k)-1
    ].

The coarse Z3 modes are k=4,8 and satisfy

    -lambda_4/rho=1

exactly.

The orthogonal modes depend on additional combinations of q_d.

Therefore the dimensionless ratios

    boxed:
    R_k
      =
    -lambda_k/rho

for k not in {0,4,8}

are observables beyond the coarse clock.

## 2. FIN values

For N=3,...,8 the reconstructed MZ shell rates give:

### N=3

    R1=1.36476
    R2=1.48525
    R3=0.96554
    R5=1.10314
    R6=1.10114.

### N=6

    R1=1.38221
    R2=1.54693
    R3=0.94307
    R5=1.17086
    R6=1.16159.

### N=8

    R1=1.41974
    R2=1.62653
    R3=0.93079
    R5=1.22047
    R6=1.20666.

These are dimensionless and remain after calibrating rho.

They also vary systematically with N.

## 3. Same-rho countermodel

To prove that these quantities are not determined by rho alone, construct for each N a positive comparator that preserves:

1. the same coarse Z3 rate k_Z3;
2. the same total 12-state exit rate;

but redistributes shell weights:
- q1,q2,q4,q5 are equalized;
- q3=q6 is chosen to preserve total exit.

This comparator has the SAME:
- rho;
- coarse Q3 law;
- total exit scale.

Yet its hidden Fourier response differs.

At the common dimensionless time

    t=1/rho,

the maximum absolute difference in hidden-mode exponential response is approximately:

    N=3:
      0.084

    N=6:
      0.103

    N=8:
      0.119.

So resolved response differs at the order of ten percent even after matching the coarse clock and total exit rate.

## 4. Experimental/operational meaning

A FIN-specific test must therefore resolve more than the Z3 coarse variable.

Candidate probes include:

- within-sector Z4/fiber relaxation;
- Fourier components of the 12 localized-state occupation pattern;
- mixed-shell response;
- pre-hydrodynamic preparation dependence;
- the short memory layer before the Z3 Markov regime.

If only the Z3 variable or long-wavelength SWAP density is measured, this information is integrated out.

## 5. Boundary

The q_d values themselves come from the conditional microscopic-to-12-state MZ reduction.

So these are FIN-specific predictions of that effective lane, not yet laboratory-confirmed constants.

But they answer an important identifiability question:

    boxed:
    FIN contains dimensionless effective observables beyond rho
    that can distinguish it from other models sharing the same coarse Q3 law.

## Verdict

P1-293 is positive.

The best near-term empirical fingerprint is not the universal diffusion law.

It is resolved pre-hydrodynamic internal response.
