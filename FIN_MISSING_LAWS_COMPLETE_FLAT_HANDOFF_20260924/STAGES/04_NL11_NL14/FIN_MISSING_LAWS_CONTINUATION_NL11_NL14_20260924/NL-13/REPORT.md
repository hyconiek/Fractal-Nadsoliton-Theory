# NL-13 — full C4 static radial wall benchmark

Status: **CONDITIONAL_NUMERICAL_FULL_C4_STATIC_WALL_FOUND_FOR_G_IDENTITY**

To test whether the radial-wall direction survives after releasing the
straight-line constraint, I declared one benchmark constitutive law:

`E[s]=int [1/2 |s_x|^2 + Phi_g_eq(s)] dx`

on the aligned C4 field.  `G=I` is **not** claimed to be sourced by FIN; NL-11
shows four D12-compatible weights remain open.

The full four-component Euler-Lagrange boundary-value problem

`s_xx = grad Phi(s)`,
`s(-infinity)=0`,
`s(+infinity)=s_localized`

was solved directly on increasingly large finite boxes.

Results:

| half-box | wall energy | lowest three linearized eigenvalues | translation correlation |
|---:|---:|---|---:|
| 15 | 0.703341373 | [0.00025131643700604335, 0.09215425838418281, 0.10299714861585522] | 0.999327 |
| 20 | 0.702879444 | [1.1823231007788633e-05, 0.08507797816398502, 0.09229402116938853] | 0.999967 |
| 25 | 0.702857096 | [-1.2473266265761748e-06, 0.08253681489209166, 0.08665765774030311] | 0.999998 |

The lowest eigenvalue tends to zero with box enlargement and its eigenvector
has >0.999 correlation with `s_x`, exactly as expected for the translation
zero mode.  The next tested eigenvalues remain positive (~0.08 and above).
The wall energy stabilizes near **0.702857096**.

The full solution is close to, but not exactly on, the simple radial line:
maximum Euclidean departure on the L=25 profile is
`0.022558`.

**Scientific meaning:** the radial coexistence discovered in NL-10 survives a
full four-component static-wall calculation under a declared positive
constitutive metric.  This is the first concrete localized field profile in
this missing-law campaign that uses the actual FIN double-well radial
potential rather than an invented sine-Gordon/DNLS potential.

It is still conditional mathematics.  The internal stiffness weights must be
sourced before this can be called a FIN-predicted wall.
