# GEOMETRY-CONVOLUTION-08

Status: EXACT RECURSION TRANSFORM; computational M=64 frontier remains open.

For a balanced split c=a+b with |c|=m and |a|=|b|=m/2, define the phase-feature
sum T(c)=sum_i c_i v_i.  The Ward split cost is exactly

Delta(a,b) = ||T(a)-T(b)||^2/m
           = [2||T(a)||^2+2||T(b)||^2-||T(c)||^2]/m.

The partition recursion uses child inverse temperature 2 beta.  Define the
Gaussian-gauged partition coefficient

Y_c(beta) = exp[- beta ||T(c)||^2/(2m)] Z_c(beta).

Then all non-diagonal Ward Boltzmann factors cancel and the exact recursion is

Y_c(beta) = exp[beta ||T(c)||^2/(2m)]/2 *
  { sum_(ordered a+b=c, |a|=m/2) Y_a(2 beta) Y_b(2 beta)
    + 1_(c even) Y_(c/2)(4 beta) }.

The second term is exactly the unordered equal-child correction.  No phase or
large-M approximation is used.

This turns the expensive split kernel into a coefficient convolution plus a
pointwise Gaussian factor.  It therefore provides an exact route to sparse
polynomial/FFT/tensor implementations for M=64 rather than further brute-force
memo recursion.

Regression of an independent C++ implementation:
- 8 phases, n=1, alpha=1: logZ = -13.888728925194343
- 8 phases, n=2, alpha=1: logZ = -13.59852846002334
- 8 phases, n=4, alpha=0.745: logZ = -6.988771497919186

These agree with the original count-vector recurrence at roughly 1e-14.

A direct n=8/M=64 memoized run was deliberately stopped after ~3 minutes rather
than promoted as a result.  At that point it had already traversed >80,000
canonical states and >139 million split terms.  M=64 therefore remains open,
but the obstruction is now computational and the convolution transform gives a
specific exact next algorithm.
