# PARAMETER-FREE-THREE-STATE-GENERATOR-142
## A memory-corrected equilateral coarse generator emerges from one microscopic process

Date: 2026-09-26

Status:
- exact finite-state slow spectrum for N=3..6;
- first-memory-moment approximation derived from the same microscopic generator;
- no fitted transition rate.

D12 symmetry makes the three localized sectors equivalent.

After the fast transition/residual mode is eliminated, the only symmetric
three-state generator has the form

    Q_eff
      =
      k *
      [ -2  1  1 ]
      [  1 -2  1 ]
      [  1  1 -2 ].

Its nonzero eigenvalues are

    -3k, -3k.

Therefore the exact slow full-process eigenvalue lambda_slow determines

    k_exact = |lambda_slow|/3.

The memory-moment closure of report 141 predicts

    k_MZ = |lambda_MZ|/3.

Results:

    N=3: k_exact=0.043813262065; k_MZ=0.044001985555; rel.err=0.4307%; N_eff(parity)=1.0038; N_eff(3-sector)=3.1293
    N=4: k_exact=0.023897425827; k_MZ=0.023951461010; rel.err=0.2261%; N_eff(parity)=1.0082; N_eff(3-sector)=3.3784
    N=5: k_exact=0.013395585305; k_MZ=0.013408252219; rel.err=0.0946%; N_eff(parity)=1.0034; N_eff(3-sector)=3.0529
    N=6: k_exact=0.007537091729; k_MZ=0.007539637385; rel.err=0.0338%; N_eff(parity)=1.0014; N_eff(3-sector)=3.1098

The error decreases from about 0.43% at N=3 to about 0.034% at N=6.

## 1. No new coarse rate parameter

The rate k is not inserted as

    c * gap

and is not calibrated independently.

It is generated from:
- the exact leave-one-out microscopic generator;
- the D12-selected basin partition;
- the zeroth and first moments of the eliminated memory kernel.

This meets an important part of the new campaign's success criterion.

## 2. Information carried by the coarse variable

At N=6 the equilibrium macro entropy gives an effective state count

    parity-basin projection:
      exp(H) ≈ 1.001432;

    three-sector projection:
      exp(H) ≈ 3.109787.

The parity projection looks dynamically easy partly because almost all
equilibrium mass is in its catch-all "other" class.

The three-sector variable actually resolves about three populated coarse states.

Thus raw semigroup error alone is a misleading model-selection criterion.

## 3. Initial slip

The first-moment low-frequency propagator contains an initial-slip amplitude.
It is not accurate at arbitrarily short times.

After the memory boundary layer, however, the slow eigenvalue is reproduced
with sub-percent accuracy.

The effective theory should therefore be stated with:
- an initial preparation/slip map;
- then the local three-state generator.

This is exactly the architecture expected from a controlled projection with
fast memory.

## 4. Current status

This is the strongest post-130 evidence so far that FIN can produce a
nontrivial effective unit from one microscopic law without adding a new
dimensionless rate parameter.

It is still a finite-N, finite-state result and not yet a recursive hierarchy.
