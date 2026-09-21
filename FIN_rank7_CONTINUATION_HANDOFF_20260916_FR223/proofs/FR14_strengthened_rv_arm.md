# FR14 — strengthened r-v arm

The FR8 anisotropic certificate can be enlarged slightly in the signed `r`
direction while keeping its other coordinates unchanged.  The common
R7P-068 interval-AD / Schur-cone checker certifies

- `|r-r_*| <= 1/6400`,
- `u <= 1/8192`,
- `v <= 1/4600`,
- `e <= 1/100000`.

Hence `lambda2(M4)<=sigma_*` throughout this box.  The same checker rejects
`|r-r_*|<=1/6300` with the other bounds unchanged; this is a method negative
control, not a physical violation.
