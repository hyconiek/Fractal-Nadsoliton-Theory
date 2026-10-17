# UE-001 source delta and nonselection

**Disposition:** `CONDITIONAL_COUPLED_FIELD_ONLY — UNVERIFIED_CURRENT_SOURCE`. The present order provides an interface and prior-scoped summary, not enough executable source data for the exact existing action.

## Necessary actual source objects

- The signed/symmetrized 12×12 strict/legacy kernel matrix `K_s` in the correct domain: **not bound**. In particular an unqualified kernel `K(d)` cannot be substituted between legacy and strict formulations.
- `M2`, `Y` and all quartic and sextic coefficient tensors with reality, symmetry, factor-of-two and potential-vs-action conventions: **not bound**.
- Diagonal/self interaction subtraction in `K_s`; kinetic signature and field dimension; boundary/Cauchy class; whether the background is a solution of **the full** action: **not bound**.
- The exact 12-complex-field plus 1-real-field action as opposed to a 13-real-coordinate mechanical polynomial restriction: **not bound**.
- Coefficient-by-coefficient strict provenance from nadsoliton: **not provided**. A moment map depending on independent `lambda` cannot fix `lambda`.

An invariant real slice is legal **only if** complex conjugation is an exact symmetry of the action and the initial field and canonical velocity/momenta lie in its fixed locus. Uniqueness of the full Cauchy evolution then retains the real slice. Merely declaring imaginary parts zero is insufficient if phase-sensitive terms or non-real couplings force them.

The fixed-kernel, same-Hessian alternatives described by the user have a general exact obstruction: if a polynomial potential `V(q)` is replaced by `V(q)+a||q-q0||^4` (with any admissible a for which the term respects the relevant stated symmetry), its value, first derivative and Hessian at `q0` agree, but its nonlinear evolution changes. This witness only establishes **underdetermination within a permitted polynomial family**; it does not license adding that term to FIN. Similarly a one-body force does not determine a two-body interaction term: a vanishing-on-single-body cross polynomial can change composed dynamics without changing isolated motion.

There is no new typed source map and no changed live dependency. **No FAR/source inventory was reopened.**
