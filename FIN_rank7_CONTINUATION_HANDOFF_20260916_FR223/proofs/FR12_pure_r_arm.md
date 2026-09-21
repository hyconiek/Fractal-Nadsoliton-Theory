# FR12 — strong pure-`r` anisotropic arm

Using the same physical local coordinates and the same interval-AD / R7P-044
Schur-cone checker as R7P-068, the following box is certified:

- `|r-r_*| <= 1/4096`,
- `u=1-exp(-3J4/2) <= 1/8192`,
- `v=1-exp(-J5/2) <= 1/8192`,
- `e=1-q_even <= 1/1024`.

Thus `lambda2(M4)<=sigma_*` throughout the box.  Relative to R7P-068 this
doubles the signed `r` radius while allowing eight times the odd-parity mass.

The same conservative checker rejects `e<=1/850` at these `r,u,v` bounds.  The
failure is a negative control for the proof method, not a physical violation.
