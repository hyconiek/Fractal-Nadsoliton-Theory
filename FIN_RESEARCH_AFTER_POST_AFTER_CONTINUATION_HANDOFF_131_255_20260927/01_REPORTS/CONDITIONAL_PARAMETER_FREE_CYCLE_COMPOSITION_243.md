# CONDITIONAL-PARAMETER-FREE-CYCLE-COMPOSITION-243
## Combining the transformation cycle with the unique swap dilation removes the continuous inter-unit coupling parameter

Date: 2026-09-26

This report combines two independently obtained structures.

## Input A: incidence candidate

Reports 183-187:

    exact deterministic information preservation
      ->
    permutation;

    one invariant orbit
      ->
    cycle C_n.

The unresolved issue is the typed interpretation of that orbit as simultaneous
physical units.

## Input B: local interaction candidate

Report 240:

    Z3 MaxEnt reset
      +
    reversible deterministic dilation
      +
    unit exchange symmetry
      +
    Z3 equivariance
      ->
    unique local gate = SWAP.

## Input C: clock/activity normalization

The isolated effective unit already supplies its refresh rate

    rho=|lambda_Z3|.

For a degree-2 cycle, requiring the same total event-participation budget per
unit and treating its two ports equivalently fixes

    r_edge=rho/2.

No independent kappa is introduced.

## Output

Conditioned on A-C, the multiunit law is fully specified:

1. vertices:
   one effective Z3 unit per cycle site;

2. topology:
   C_n;

3. local event:
   swap endpoint Z3 states;

4. edge rate:
   rho/2;

5. stationary measure:
   product uniform, decomposed by global count sectors;

6. density spectrum:

       lambda_m
         =
       rho[
         1-cos(2 pi m/n)
       ];

7. long scale:

       lambda_1
         ~
       2 pi^2 rho/n^2.

So the previous free continuous inter-cell coupling has disappeared.

## What remains genuinely unsourced

This is NOT yet a FIN theorem because three role-transfer premises remain.

### 1. Simultaneity bridge

Why are vertices of the information-preserving transformation orbit
simultaneous units rather than successive history states?

### 2. One-orbit source

Why does the multicell transformation have exactly one orbit?

### 3. Activity-budget principle

Why must composition preserve the isolated per-unit refresh activity rather
than add a new interaction activity?

If these three premises are derived, the 1D cycle composition would be
parameter-free at the dimensionless level.

This sharply narrows the old composition problem.
