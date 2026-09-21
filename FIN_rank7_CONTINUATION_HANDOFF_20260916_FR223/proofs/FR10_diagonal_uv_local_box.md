# FR10 — diagonal `u-v` anisotropic local certificate

This continuation keeps the completed R7P ledger unchanged and enlarges only the
locally certified physical four-amplitude neighborhood around the R7P-068 equality
point.

Use the same local coordinates as R7P-068:

- `x=r-r_*`, `r=exp(-2 J3)`, signed;
- `u=1-exp(-3 J4/2)`;
- `v=1-exp(-J5/2)`;
- `e=1-q_even`.

The original certificate used the common radius `1/8192` in all four coordinates.
FR10 runs the identical second-order rational interval-AD and Schur-cone proof on
an anisotropic box.

## Certified box

The checker certifies

`|x| <= 1/8192`, `0<=u<=1/4800`, `0<=v<=1/4800`, `0<=e<=1/3072`.

Throughout this box the six sign prerequisites remain strict, `c2>0`, and both
Schur endpoint cone tests pass after the same `1/10000` margin used by R7P-068.
Therefore the R7P-044 disjunction

`P<=0 OR P1>=0`

holds throughout the box, and hence

`lambda2(Mtilde)<=sigma_*`, equivalently `lambda2(M4)<=sigma_*`.

This is a simultaneous enlargement in both off-face coordinates `u,v` and in
odd-parity mass `e`; it is not obtained by taking a union of one-coordinate arms.

## Negative controls

The same conservative checker rejects both

- `e<=1/3052` with `u,v<=1/4800`, and
- `u,v<=1/4700` with `e<=1/3072`.

These failures are deliberately retained as checker-conservatism controls. They
are **not** counterexamples to the ceiling.

The theorem remains local to the physical shared-field four-amplitude model and
does not imply a full positive-orthant theorem or any full-seven-coordinate
Hessian statement.
