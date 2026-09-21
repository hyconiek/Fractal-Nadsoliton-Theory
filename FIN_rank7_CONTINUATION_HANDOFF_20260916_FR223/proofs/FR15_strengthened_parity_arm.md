# FR15 — strengthened FR9 parity arm

A fixed-seed search on the compact residual repeatedly approached the artificial intersection of the old FR9 `e=1/100000` wall and the FR11 `u=1/4800` wall, always with negative physical gap.  Replaying the same second-order interval-AD / R7P-044 checker on the full FR9 geometry gives a much larger safe parity radius:

- `|r-r_*| <= 1/6500`,
- `u <= 1/2432`,
- `v <= 1/8192`,
- `e=1-q_even <= 1/14000`.

The checker certifies the full box.  As a proof-method negative control, the otherwise identical box with `e<=1/13000` fails the boundary Schur test.  This failure is not a physical counterexample.

Thus FR15 strictly subsumes FR9 and removes the residual-search wall that motivated it.
