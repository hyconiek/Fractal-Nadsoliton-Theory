# T52 Current Ideal Isotropic Void Limit Omega-Phi Theta-Reduction Nonadmissibility Boundary Spec

Status: `T52_CURRENT_IDEAL_ISOTROPIC_VOID_LIMIT_OMEGA_PHI_THETA_REDUCTION_NONADMISSIBILITY_BOUNDARY_SPEC_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `N316`, one stronger proposal appears:

```text
force the candidate route into the idealized special case
  A_pair^cand = I_2
  lambda_1 = lambda_2 = 1
and then read off
  theta_1 = arctan(omega)
  theta_2 = arctan(phi)
as an actual theta reduction
```

The question is not whether this specialization is easy to write.

It is easy to write.

The narrower question is:

```text
can the current repo honestly promote that idealized specialization
from candidate-law notation
to actual theta reduction
without false pass?
```

## Scope

`T52` is scoped only to this proposed idealized limit on the current
omega-phi primordial-preorientation route.

It reuses only:

1. `N298`
2. sandbox `N18`
3. `N314`
4. `N315`
5. `N316`
6. `R1/C48/C49`
7. `H42/C29` only as a check that unrelated appearances of `I_2` do not
   discharge this route.

It does **not** decide:

1. impossibility in principle of future isotropic limit routes,
2. impossibility in principle of future actual theta reduction,
3. impossibility in principle of component 2,
4. impossibility in principle of future loop breaking.

## Exact proposal under audit

The proposed specialization is:

```math
A^{cand}_{pair}=I_2
```

```math
\lambda_1=\lambda_2=1
```

so that the candidate map-rule from `N316` collapses formally to:

```math
\boldsymbol{\vartheta}^{cand}_{pair}
=
\mathcal{N}_{pair}
\begin{pmatrix}
\omega\\
\phi
\end{pmatrix}
```

and, under one explicit normalization choice,

```math
\theta_1 = \arctan(\omega),
\qquad
\theta_2 = \arctan(\phi).
```

## Exact decision question

Can the current repo honestly export anything stronger than:

```text
one actual packaged pair-map-rule candidate exists
```

namely:

```text
one actual theta-reduction result obtained by forcing
the ideal isotropic void limit
  A_pair^cand = I_2
  lambda_1 = lambda_2 = 1
```

on the current repo state?

## Boundary target

If the answer is negative, freeze:

```text
IdealIsotropicVoidLimit_omega_phi_theta_actual_reduction_nonadmissibility_boundary_v1
```

with the intended meaning:

```text
the ideal isotropic void limit may be writable as a future idealization,
but the current repo exports no route-local discharge
for A_pair^cand = I_2,
no route-local discharge for lambda_1 = lambda_2 = 1,
no actual theta values,
and no actual theta reduction;
therefore that idealized specialization is not currently admissible
as an actual theta-reduction step
```

## Hard limits

`T52` must not claim:

1. actual `A_pair = I_2`,
2. actual `lambda_1 = lambda_2 = 1`,
3. actual `theta_1 = arctan(omega)`,
4. actual `theta_2 = arctan(phi)`,
5. actual theta reduction,
6. actual pair population,
7. actual component-2 support,
8. actual sandbox loop break,
9. actual `E_orient`,
10. admissible `S_sel_int`,
11. strict-core selector closure,
12. `QW-2191` discharge,
13. ToE closure.
