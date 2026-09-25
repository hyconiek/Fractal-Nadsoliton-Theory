# TREE-DYN-INVERSE-05 — reconstructing the dynamic tree

Status: **CONSTRUCTIVE_CONDITIONAL_THEOREM**

## Main result
The static boundary response gives an effective-resistance additive tree metric, which identifies the minimal positive weighted tree up to suppressed degree-2 vertices.  The high-frequency expansion of the dynamic response then exposes leaf-parent storage coefficients; stripping terminal layers recursively reconstructs internal storage in the declared class.  Pole positions alone are insufficient; matrix residues/high-frequency coefficients are necessary.

## Core formulas
\[
d_{ij}=(e_i-e_j)^T\Lambda(0)^+(e_i-e_j),
\]
\[
\Lambda(z)=L_{BB}-z^{-1}L_{BI}C_I^{-1}L_{IB}+O(z^{-2}),
\]
\[
c_v=rac{g_i^2}{-\lim_{z	o\infty}z[\Lambda(z)-L_{BB}]_{ii}}.
\]

## Evidence / reproduction
A two-internal-node counterexample with identical generalized eigenvalues but different storages demonstrates that poles without residues are non-identifying.

## Caveats
Dynamic degree-2 subdivisions and nonminimal realizations require separate equivalence-class analysis.

## Next question
Derive finite-noise stability bounds and classify degree-2 dynamic ambiguities.
