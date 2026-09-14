# R7P-016 — asymptotic/domain chart rules

Positive four-amplitude fields use the exact shared field
`h_j=J3 cos(alpha_j)+J4 cos(beta_j)+J5 cos(alpha_j-beta_j)+J6 cos(2 alpha_j)`.
A large-field limit is controlled by exposed supports of the 12 C4 feature
points with a normal vector in the nonnegative orthant. The receiving code
records candidate exposed supports numerically, but they are not promoted to
exact boundary strata until certified.

For proof charts, use relative softmax gaps `delta_j=max_i h_i-h_j>=0` and
weights `z_j=exp(-delta_j)`. This preserves all shared-field relations.

Do **not** polynomialize the model by treating `exp(J5/2)` and
`exp(sqrt(3) J5/2)` as independent variables. Odd-sector mode-5 values contain
`sqrt(3)/2`; independent exponentials would enlarge the physical domain.
Where such irrational powers occur, retain a validated transcendental interval
chart or an exact algebraic relation with an explicit shared field variable.

These positive-amplitude charts do not cover arbitrary signed theta7 limits;
full-seven-coordinate asymptotics are a separate domain problem.
