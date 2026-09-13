# R7P-034 — invariant two-harmonic stationary system

Take

`h_j = J cos(pi j/2) + K (-1)^j`.

The field has period four and is reflection-even. Hence its Gibbs probability
has period four, and all Fourier expectations outside sectors `0,3,6,9` vanish;
the sine component of sector 3 also vanishes. Thus this is an exact invariant
stationary family of the full X7 equations, not merely a C4 truncation.

Summing the four residue classes gives

`x = E[cos(pi j/2)] = sinh(J)/(cosh(J)+exp(-2K))`,

`y = E[(-1)^j] = (cosh(J)-exp(-2K))/(cosh(J)+exp(-2K))`.

With the normalized X7 convention, full stationarity is therefore exactly

`J = g (lambda3/6) x`,
`K = g (lambda6/12) y`.

The compact search rectangle follows without approximation from `|x|,|y|<=1`:
`|J|<=g lambda3/6`, `|K|<=g lambda6/12`.
