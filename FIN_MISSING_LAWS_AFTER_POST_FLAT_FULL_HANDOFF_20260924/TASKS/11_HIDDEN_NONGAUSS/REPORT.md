# First non-Gaussian hidden-memory correction

Status: **CONDITIONAL_KRAMERS_MOYAL_RESULT**

## Result
For x=sqrt(N)(mu-mubar), the third Kramers-Moyal tensor enters at N^-1/2.
It is affine in the four hidden amplitudes, so substituting the exact hidden
solution produces an explicit exponential history dependence. At equilibrium
the third tensor itself vanishes by j<->k symmetry, but its hidden derivatives
do not.

## Key formulas
\[C^{(3)}=C_0^{(3)}(\mu)+\sum_{a=1}^4\nu_aT_a^{(3)}(\mu),\quad
\nu_a(t)=e^{-t}\nu_a(0)+\int_0^te^{-(t-s)}h_a(\mu_s)ds.\]

## Caveat
This is the first non-Gaussian correction, not a physical memory time in SI units.

## Next question
Derive the full first Edgeworth operator after hidden elimination.
