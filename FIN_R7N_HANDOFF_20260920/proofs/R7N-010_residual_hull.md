# R7N-010 — residual hull and tails

On the nonnegative-field chart
`r=exp(-2J3), s=exp(-3J4/2), t=exp(-J5/2), y=exp(-2J6)`, every finite point lies in `(0,1]^4`.
The retained compact hull is
`r>=1/900, s>=1/128, t>=1/9, y>=1/1000000`.

Outside it, at least one strict lower inequality holds. If `r<1/900`, then
`exp(-J3)=sqrt(r)<1/30`; if `s<1/128`, the accepted FR1 `J4` tail applies;
if `t<1/9`, the accepted FR1 `J5` tail applies; if `y<1/1000000`, the accepted
FR42 large-`J6` tail applies. The lower endpoints themselves are intentionally
retained in the compact hull, so tail/hull overlap is allowed and no boundary is removed.
The infinite-field faces are reached by continuity of the normalized weight law,
whose anchor weight is exactly one.

This is a domain-decomposition lemma conditional on the accepted FR1/FR42 tail
certificates. It does **not** assert the remaining hull is safe.
