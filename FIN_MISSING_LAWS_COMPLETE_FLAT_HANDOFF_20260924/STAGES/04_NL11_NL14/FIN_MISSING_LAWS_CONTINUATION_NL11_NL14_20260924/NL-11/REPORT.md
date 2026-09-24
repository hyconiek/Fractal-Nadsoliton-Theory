# NL-11 — radial intercell stiffness source

Status: **PARTIAL_SOURCE_D12_REDUCES_RADIAL_STIFFNESS_TO_FOUR_WEIGHTS_BUT_DOES_NOT_SELECT_THEM**

The internal aligned C4 amplitudes are not four copies of one scalar.  In the
full X7 representation, `k=3,4,5` are inequivalent two-dimensional D12
irreducible sectors and `k=6` is the one-dimensional parity sector.

Therefore a D12-invariant positive quadratic form has the block structure

`G = diag(g3 I2, g4 I2, g5 I2, g6)`, `gk>0`.

On the reflection-even C4 chart this becomes simply

`diag(g3,g4,g5,g6)`.

So symmetry removes all cross-sector mixing and reduces a generic symmetric
4x4 tensor from 10 parameters to **four positive weights**.

The NL-02 spatial law can then be extended without changing its refinement
algebra:

`E_ij = 1/2 c_ij (theta_i-theta_j)^T G (theta_i-theta_j)`,

with `c_ij=kappa Area(F_ij)/ell_ij` (and `kappa/ell` on S1).

But D12 does not choose the four weights.  Three immediately available
D12-invariant choices are already inequivalent:

- identity weights;
- feature-Gram weights `(lambda3,lambda4,lambda5,lambda6)`;
- A-weighted weights `(lambda3^2,...,lambda6^2)`.

For the current strict operator the four lambdas are
`1.961406862, 2.199568849, 2.298606272, 2.342182041`.

**Conclusion:** geometry and internal constitutive response separate.  FIN now
reduces the missing radial stiffness law to four sector coefficients, but does
not yet source them.
