# MP7-014 — no nonzero boundary branch for `0<g<=250/67`

The only structurally admissible proper nonempty supports from MP7-013 are `{6}`, `{4}`, `{4,6}`, `{3,6}`.  All are excluded throughout the Target-P gain window.

Write `a3^2=lambda3/6`, `a4^2=lambda4/6`, `a6^2=lambda6/12` and fields `Jk=a_k s_k`.

## Pure `{6}`

The stationary equation is

`J6 = g (lambda6/12) tanh(J6)`.

For `J6>0`, `tanh J6 < J6`, while the strict spectral upper endpoint gives

`g lambda6/12 <= (250/67) lambda6_hi/12 < 1`.

Thus no nonzero solution exists.

## Pure `{4}`

With `x=3J4/2`,

`E[cos(2*pi*j/3)] = (e^x-1)/(e^x+2) =: m(J4)`.

For all `J4>=0`,

`m(J4) <= (2/3) J4`.

Indeed the inequality is equivalent to

`F(x)=4x(e^x+2)-9(e^x-1)>=0`.

Here `F(0)=0`, `F'(x)=e^x(4x-5)+8`, and `F''(x)=e^x(4x-1)`.  The derivative has its minimum at `x=1/4`; the elementary series bound `e^(1/4)<2` gives `F'(1/4)>0`, hence `F` is increasing and nonnegative.

Stationarity would give

`J4 = g(lambda4/6)m(J4) <= g lambda4 J4/9`,

but `(250/67)lambda4_hi/9 <1`, contradiction for `J4>0`.

## `{4,6}`

Modes four (period three) and six (parity) are independent under the uniform `Z12` label average because `Z6 ~= Z3 x Z2`.  Explicitly the partition factorizes into the pure-mode-4 factor times `cosh J6`.  Consequently the two mean equations decouple into the already excluded pure `{4}` and `{6}` equations.  No nonzero `{4,6}` stationary branch exists in the window.

## `{3,6}`

Over `j mod4`,

`Z = [e^{J6} cosh J3 + e^{-J6}]/2`,

so

`E[C3_unscaled] = sinh J3/(cosh J3+e^{-2J6})`,
`E[C6_unscaled] = (cosh J3-e^{-2J6})/(cosh J3+e^{-2J6})`.

The second equation implies

`0<=J6<=g lambda6/12`, hence

`2J6 <= (250/67)lambda6_hi/6 < 3/2`.

A rational exponential-series bound gives `e^(3/2)<9/2`, so

`c=e^{-2J6} > 2/9`.

For `x>=0` and `0<c<=1`,

`sinh x/[x(cosh x+c)] <= 1/(1+c)`.

To see this, the desired inequality is
`x cosh x-sinh x >= c(sinh x-x)`; it is enough to use `c<=1` and the stronger
`x cosh x-sinh x >= 2(sinh x-x)`, equivalent to
`x cosh x-3 sinh x+2x>=0`.  Three derivatives reduce this to `x sinh x>=0` with zero initial data.

Therefore for `J3>0`,

`E[C3_unscaled]/J3 < 9/11`.

The mode-three stationary equation would require

`1 = g(lambda3/6) E[C3_unscaled]/J3
   < (250/67) lambda3_hi * 3/22 < 1`,

a contradiction.  The final strict rational margin is recorded in `results/MP7-014_boundary_gain_exact_checks.json`.

## Conclusion

For every `0<g<=250/67`, the only aligned nonnegative stationary support on the boundary is the uniform zero-amplitude point. Any nonzero aligned stationary point in this gain window must have **all four C4 amplitudes strictly positive**.

This does not assert uniqueness of the interior branch or global-minimum uniqueness.
