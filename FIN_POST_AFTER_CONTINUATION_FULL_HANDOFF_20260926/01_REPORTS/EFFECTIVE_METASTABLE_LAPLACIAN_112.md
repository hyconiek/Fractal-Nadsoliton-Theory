# EFFECTIVE-METASTABLE-LAPLACIAN-112
## D12 symmetry gives an explicit circulant transition operator for the twelve-minimum reduction

Date: 2026-09-26

Status:
- exact representation formula once one assigns one symmetric rate k_d to each
  pair-saddle class;
- only exponential rate exponents are fixed by the current FIN results;
- prefactors remain open.

Let k_d be the effective transition rate per edge between minima at cyclic
separation d.

D12 symmetry makes k_d independent of the starting label.

For d=1,...,5 there are two neighbors ±d.
For d=6 there is one antipodal neighbor.

The effective generator on functions f:Z12->R is

    (Qf)(n)
      =
      sum_(d=1)^5
        k_d[
          f(n+d)+f(n-d)-2f(n)
        ]
      +
        k_6[
          f(n+6)-f(n)
        ].

The Fourier modes

    exp(2 pi i m n/12)

diagonalize Q.

The positive decay rates are

    Lambda_m
      =
      2 sum_(d=1)^5
        k_d[1-cos(2 pi m d/12)]
      +
      k_6[1-(-1)^m].

At large g,N,

    -log k_d / N
      =
      alpha_d g-log2+o(g)

at the exponential level.

Since alpha_3 is uniquely smallest,

    k_3 >> k_4 >> all remaining k_d

exponentially.

The d=3 operator alone has a three-dimensional zero space:
- the global constant;
- two modes distinguishing the three residue classes modulo 3.

The next rate k_4 lifts precisely those inter-sector modes.

For m=4 and m=8,

    Lambda_m
      =
      3 k_4
      + exponentially smaller corrections

in the hierarchy k3>>k4>>....

Thus, provided the subexponential prefactor of k4 is nonzero,

    spectral-gap exponent
      =
      alpha_4 g-log2.

By contrast, relaxation inside each four-state d3 cycle occurs on the faster
k3 scale.

Hence the effective dynamics has two separated relaxation stages:

    fast:
      equilibration inside each 4-state sector;

    slow:
      equilibration among the three sectors.

This is a concrete dynamical consequence of the stationary support/barrier
classification.

What is still missing for a physical prediction:
- the Eyring-Kramers/potential-theory prefactors;
- the conversion from FIN refresh time to seconds;
- evidence that these twelve minima correspond to physical configurations.
