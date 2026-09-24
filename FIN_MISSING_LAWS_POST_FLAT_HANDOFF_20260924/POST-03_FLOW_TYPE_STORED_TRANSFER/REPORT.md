# POST-03 — FLOW-TYPE / stored transfer

Status: **EXACT_MINIMAL_LIFT_PLUS_SOURCE_GAP**

For the moment-conditioned edge covariance `D_c` and visible covariance
`A7=B D_c B^T`, the minimum-seminorm edge realization of a visible mediator
state `h` is

`Q*(h)=D_c B^T A7^+ h`.

Replay:
- `||B Q* - h|| = 1.459e-14`;
- edge minimal cost and visible `h^T A7^+ h` differ by
  `7.105e-15`.

Thus A7^+ is the quotient metric induced by the parent edge-flow geometry.

But Q* is fully determined by h.  A genuinely larger state is

`Q = Q*(h) + z`, `z in ker B`.

For complete C12, `dim ker B = 55`.  Only z can store history that is not
already contained in h.

If Q is *declared* to be accumulated transfer and h=BQ, then with J=-Qdot the
continuity equation `hdot + BJ=0` is a kinematic identity.  FIN still does not
source why Q rather than only h is physical, nor an evolution law for z.
