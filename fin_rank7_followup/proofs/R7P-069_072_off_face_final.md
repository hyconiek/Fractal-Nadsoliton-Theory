# R7P-069--072 — terminal partial off-face result

The campaign does **not** prove the global positive-orthant four-amplitude
ceiling.  It does prove several exact/certified subdomains and leaves one
explicit residual compact core.

## New positive-volume tail

Put `t=exp(-J5/2)`.  Because all four fields are nonnegative, relative to label
`j=0` every non-anchor Boltzmann weight is bounded above by
`t^(2(1-cos(5 pi j/6)))`; the other nonnegative fields can only decrease that
ratio.  For the four C4 feature vectors `F_j`,

`tr Cov(F) <= E ||F-F_0||^2 <= sum_{j!=0} ratio_j ||F_j-F_0||^2`.

Using the accepted outward strict spectral intervals and interval arithmetic at
`t=2^-11`, the complete sum is `<0.5050012544110`, whereas
`2 sigma_*>0.5348864884576`.  Hence throughout the entire tail
`t<=2^-11`, independently of `J3,J4,J6`,

`lambda2(M4) <= tr(M4)/2 < sigma_*`.

The replay is `PYTHONPATH=.:..:src python src/off_face_tail.py`.

## Other certified pieces

The exact boundary-Ising closure is covered by R7P-055; the shared-field
intraparity covariance by R7P-063; the extreme face has the exact face theorem;
and R7P-068 supplies a finite local physical-cone box of radius `1/8192` around
the compactified equality point.

## Residual and strongest theorem

The remaining compactified core with `t>2^-11`, after removing those certified
faces/boundaries/local boxes, has not received a complete verified cover.  The
bounded adversarial search R7P-066 found no admissible violation but is not a
proof.  Therefore R7P-070 records **no certified counterexample**, R7P-071
exports only the union-of-certified-subdomains theorem, and R7P-072 forbids
promotion to a full positive-orthant 4D ceiling or to the already-refuted
full-seven-coordinate everywhere index bound.
