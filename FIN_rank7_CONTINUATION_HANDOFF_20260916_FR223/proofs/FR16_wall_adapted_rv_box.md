# FR16 — wall-adapted `r-v` box

After FR15, every fixed-seed residual search moved to the old FR14 `|r-r_*|=1/6400` wall at `v≈1/5654`, with a negative physical gap near `-6.13e-6`.

Instead of extending the full FR14 box uniformly, narrow the `v` extent to match that wall.  The same interval-AD / characteristic-inertia checker certifies

- `|r-r_*| <= 1/5000`,
- `u <= 1/8192`,
- `v <= 1/5600`,
- `e <= 1/100000`.

This contains the numerical wall point and increases the certified signed `r` radius by about 28% relative to FR14.  The otherwise identical box with `|r-r_*|<=1/4800` fails the conservative boundary Schur test; that is a method negative control, not a physical violation.
