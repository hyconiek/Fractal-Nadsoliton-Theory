# EXACT-Z3-LUMPING-OF-LOCALIZED-GENERATOR-217
## The 12-state circulant generator lumps exactly to the previously derived three-state coarse law

Date: 2026-09-26

Status:
exact algebraic consequence of circulant shell rates, numerically instantiated
for N=3..8.

Partition the localized labels by

    j mod 3.

From any state in one sector, the one-direction rate to either of the other
Z3 sectors is

    boxed:
    k_Z3
      =
      q1+q2+q4+q5.

The pure fiber/within-sector shells are

    d3,
    d6.

Therefore the exact lumped three-state generator is

    Q_Z3
      =
      k_Z3
      [ -2 1 1
         1 -2 1
         1 1 -2 ].

Its nonzero eigenvalue is

    lambda_Z3=-3 k_Z3,

which agrees identically with the k=4,8 Fourier eigenvalue of the 12-state
generator.

Results:


    N=3:
      k_Z3=0.0440019855553
      direct d4 fraction=48.023 %
      d1 fraction=4.069 %
      d2 fraction=21.182 %
      d5 fraction=26.726 %

    N=4:
      k_Z3=0.0239514610098
      direct d4 fraction=48.776 %
      d1 fraction=5.299 %
      d2 fraction=20.009 %
      d5 fraction=25.915 %

    N=5:
      k_Z3=0.013408252219
      direct d4 fraction=49.884 %
      d1 fraction=5.961 %
      d2 fraction=18.958 %
      d5 fraction=25.198 %

    N=6:
      k_Z3=0.00753963738482
      direct d4 fraction=51.307 %
      d1 fraction=6.215 %
      d2 fraction=17.960 %
      d5 fraction=24.518 %

    N=7:
      k_Z3=0.00421397035429
      direct d4 fraction=53.033 %
      d1 fraction=6.144 %
      d2 fraction=16.987 %
      d5 fraction=23.837 %

    N=8:
      k_Z3=0.00232997400889
      direct d4 fraction=54.963 %
      d1 fraction=5.874 %
      d2 fraction=16.031 %
      d5 fraction=23.132 %

## Important correction to a single-saddle picture

At N=8 the pure d4/base jump supplies about 55% of the one-direction
inter-sector rate.

The other ~45% comes from mixed shells d1,d2,d5.

Thus the successful Z3 coarse rate is not yet controlled by the d4 saddle
alone.

The quotient is exact because of symmetry, while its effective rate is a sum
over several microscopic/mesoscopic transition channels.

This is precisely why a capacity/committor treatment is preferable to assigning
one rate from one barrier by inspection.

## Effective capacity

For the symmetric three-state generator:

    pi(C0)=1/3,

and the total exit rate is

    2 k_Z3.

Hence

    cap_eff(C0,C1 union C2)/pi(C0)
      =
      2 k_Z3
      =
      (2/3)|lambda_Z3|.

This is exactly the long-time rate used in reports 142/145.

So the memory-renormalized generator resolves the earlier mismatch:
raw microscopic whole-basin capacity contains recrossing flux,
whereas capacity of the effective metastable chain equals the committed
long-time transition rate by construction.
