# GEOMETRY-TRANSITION-09 — exact-convolution M=64 extension

Status: **RECOMPUTED NUMERICAL FINITE-SIZE RESULT; NO THERMODYNAMIC THEOREM**.

The exact Gaussian-gauged count-vector recursion from GEOMETRY-CONVOLUTION-08
was evaluated for the symmetric eight-phase alphabet at M=64 (n=8 replicas of
each phase).  The convolution at the M=32 layer is performed by FFT, but the
mathematical recursion is unchanged.  Selected convolution coefficients were
checked against direct summation to ~2e-15, and every dense M=64 point used a
same-alpha M=32 regression against the independent direct C++ recurrence.

Using C_H = alpha^2 d^2(log Z)/d alpha^2, the local maxima are:

- M=16 (n=2): alpha*=0.7341336, C_H,max=11.09167;
- M=32 (n=4): alpha*=0.7459035, C_H,max=63.17082;
- M=64 (n=8): alpha*=0.82457, C_H,max=300.607.

The M=64 value comes from a 0.0005 alpha grid around the maximum; a local
quadratic interpolation gives alpha*=0.8245696 and C_H,max=300.60694.

The peak-height doubling ratios are about 5.70 and 4.76 (local exponents 2.51
and 2.25 with respect to n).  They are not stable enough to promote an
asymptotic exponent.  More importantly, the peak location itself drifts strongly
between M=32 and M=64, so the earlier impression of alpha*=~0.74 stabilization
is superseded.

Central energy cumulants obtained from derivatives of log Z at the peak give:

- M=16: Binder-like U4 ~= 0.1733, standardized skewness ~=0.6005;
- M=32: U4 ~= 0.2884, skewness ~=0.2517;
- M=64: U4 ~= 0.3923, skewness ~=0.1136.

For M=64 the fourth-derivative result is stable between h=0.0005 and h=0.001
(U4 0.39236 versus 0.39226).  The trend is compatible with an increasingly
symmetric broad/two-sector competition, but three sizes do not establish a
first-order thermodynamic transition.

At alpha=0.8245 the ordered root-split distribution remains highly nontrivial:
its Shannon effective count is about 1.05e5 splits, the largest single split
has probability <0.01, and the strongly segregated root sector q^2>=96 has
probability about 0.368.  Thus the peak is not a collapse onto one hierarchy.

The replicated-minimum theorem gives E_min(n)=n E_min(1), but minimum
degeneracy has not been proved for n=8.  Using only the four explicit symmetry
copies gives a lower bound P_min >= about 0.0064 at the M=64 peak, not an exact
minimum probability.

## Exact root-order profile (correction to an intermediate interpretation)

A local-maximum finder initially omitted the endpoint q^2=q^2_max and therefore
made the M=64 profile look approximately two-sector.  The full occupied-support
histogram corrects that interpretation.

At the M=64 peak (alpha=0.8245 used for the histogram), the endpoint q^2=128
itself carries about 0.2153 total ordered-root probability and is the largest
macrobin.  Additional local maxima occur at q^2=96 (~0.0985), 104 (~0.0350),
88 (~0.0390), 72 (~0.0413), 64 (~0.0420), etc.  The correct description is a
multimodal landscape with a strong segregated endpoint sector, not a clean
binary phase mixture.

For a reproducible coarse barrier diagnostic, compare the local q^2=96 sector
to the endpoint q^2=128 along occupied q^2 bins.  The smallest occupied-bin
probability between them is ~7.789e-4 at q^2=100, giving log probability
barriers about 4.84 from the q^2=96 side and 5.62 from the endpoint side.
At M=32 the analogous q^2=24 versus endpoint q^2=32 pair has its intervening
minimum at q^2=26, with barriers only about 1.08 and 1.37.  M=16 has no such
intervening occupied-support valley: the root-order histogram rises toward its
endpoint.  This is finite-size evidence for sharpening sector separation, not a
thermodynamic free-energy barrier theorem.
