# P297 Current Actual Ideal Isotropic Void Limit Omega-Phi Theta-Reduction Nonadmissibility Boundary Probe

Status: `P297_EXECUTED_CURRENT_ACTUAL_IDEAL_ISOTROPIC_VOID_LIMIT_OMEGA_PHI_THETA_REDUCTION_NONADMISSIBILITY_BOUNDARY_PROBE_NO_FALSE_PASS`
As of: `2026-03-09`

## Probe question

Can the current repo honestly upgrade the `N316` candidate law into actual
theta reduction by forcing the idealized special case:

```math
A^{cand}_{pair}=I_2,
\qquad
\lambda_1=\lambda_2=1
```

and then reading:

```math
\theta_1=\arctan(\omega),
\qquad
\theta_2=\arctan(\phi)
```

as an actual exported result?

## Probe checks

### Check 1: does the route already export actual values for `lambda_1`,
`lambda_2`?

NO.

`F203/N314` explicitly keep:

1. `actual_lambda_1_value_present = false`,
2. `actual_lambda_2_value_present = false`.

So the proposal cannot inherit `lambda_1 = lambda_2 = 1` from current route
data.

### Check 2: does the route already export actual discharge of
`A_pair = I_2`?

NO.

`F205/N316` export only:

```text
A_pair_candidate_in_GL2_plus = true
```

That is strictly weaker than:

```text
A_pair = I_2
```

and no route-local equality discharge is exported.

### Check 3: do unrelated occurrences of `I_2` elsewhere discharge this
route-local equality?

NO.

`H42/C29` contain unrelated isotropic or local-projector appearances of
`I_2`, but those are not:

1. the present omega-phi pair-map route,
2. a route-local `A_pair` identity witness,
3. an actual theta-reduction export.

### Check 4: does the current repo already export actual theta values on this
route?

NO.

`N316` remains explicit:

1. `actual_theta_reduction` is absent,
2. `actual_theta_export` is absent,
3. `actual_pair_population` is absent.

### Check 5: would the proposed specialization overturn `N298`?

NO.

`N298` remains correct because the current repo still exports no actual typed
map `(omega,phi) -> (theta_1,theta_2)`.

A writable idealization is not the same thing as an actual exported map.

### Check 6: would the proposed specialization overturn sandbox `N18`?

NO.

Sandbox `N18` remains correct because:

1. no actual theta values are exported by this proposal,
2. no actual populated basis-pair instance follows from it,
3. no noncyclic anchor is actually discharged.

## Probe verdict

The strongest honest result is:

```text
the ideal isotropic void limit is writable only as a symbolic special case
of the current candidate law,
but is not admissible on current repo state
as actual theta reduction
```

## Product

Export:

```text
IdealIsotropicVoidLimit_omega_phi_theta_actual_reduction_nonadmissibility_boundary_v1
```

and keep explicit:

1. `ideal_specialization_writable_only = true`,
2. `route_local_actual_identity_discharge_absent = true`,
3. `route_local_actual_lambda_values_absent = true`,
4. `actual_theta_values_absent = true`,
5. `actual_theta_reduction_absent = true`,
6. `actual_loop_break_absent = true`.
