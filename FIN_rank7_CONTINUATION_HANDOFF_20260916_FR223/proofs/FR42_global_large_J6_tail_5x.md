# FR42 — fivefold upgrade of the global large-`J6` tail

FR5 certified `y=exp(-2J6)<=1/5,000,000` by combining a boundary reserve
`delta=10^-6` with the physical odd-mass bound and the global feature-diameter
perturbation estimate.  FR42 repeats the same rigorous argument at the larger
boundary shift

`delta = 1/200000 = 5*10^-6`.

The boundary Bernstein cover now quarantines any difficult leaf using the
already-certified **FR9** or **FR16** projection.  It terminates with 637
leaves and no unresolved cells:

- 440 `SAFE_A_NONPOS`,
- 150 `SAFE_B_NONNEG`,
- 37 `LOCAL_FR9`,
- 10 `LOCAL_FR16`.

For `y<=1/1,000,000`, the exact physical parity relation gives

`e <= y/(1+y) <= 1/1,000,001`.

The FR5 diameter estimate `D^2<5` therefore gives

`||M4-C_+|| <= 5e <= 5/1,000,001 < 5*10^-6 = delta`.

The local leaves are also safe at finite `y`, because both FR9 and FR16 allow
`e<=1/100000`, much larger than `1/1,000,001`.

Hence the full physical four-amplitude tail

`exp(-2J6)<=1/1,000,000`

satisfies `lambda2(M4)<=sigma_*`.  This is a **factor-five** enlargement of the
FR5 global tail.
