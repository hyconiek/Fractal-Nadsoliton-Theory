# P1974 / S924 — Strict Bianchi-I Anisotropic EOM Residual Obstruction Theorem Packet

Status: `OPEN_OBSTRUCTION_WITH_TRACE`  
Route: `strict_only`  
Date anchor: `2026-05-18`

## Professor-level decision

`P1973` proved a local scalar finite-part transport residual between FRW and
Bianchi-I.  That was necessary but not sufficient.  The next honest scientific
question is whether the **tensor metric EOM components** also transport by that
scalar lock.  `P1974` answers: generically no.

This is a negative theorem packet.  It prevents the previous scalar transport
success from being over-promoted into background-independence closure.

## Setup

Use diagonal Bianchi-I directional Hubble rates

```text
H_1 = H + sigma_1
H_2 = H + sigma_2
H_3 = H + sigma_3
```

with trace-free anisotropy

```text
sigma_1 + sigma_2 + sigma_3 = 0,
dsigma_1 + dsigma_2 + dsigma_3 = 0.
```

Thus:

```text
sigma_3  = -sigma_1 - sigma_2
dsigma_3 = -dsigma_1 - dsigma_2.
```

## Einstein component comparison

For diagonal Bianchi-I:

```text
G_00 = H_1 H_2 + H_1 H_3 + H_2 H_3,
G_ii/a_i^2 = -(dot H_j + dot H_k + H_j^2 + H_k^2 + H_j H_k),  j,k != i.
```

For FRW:

```text
G_00^FRW = 3 H^2,
G_ii^FRW/a^2 = -(2 dot H + 3 H^2).
```

Subtracting FRW from Bianchi-I gives the exported residual vector:

```text
R_00 = -sigma_1^2 - sigma_1 sigma_2 - sigma_2^2
R_11 =  3H sigma_1 + dsigma_1 - sigma_1^2 - sigma_1 sigma_2 - sigma_2^2
R_22 =  3H sigma_2 + dsigma_2 - sigma_1^2 - sigma_1 sigma_2 - sigma_2^2
R_33 = -3H sigma_1 - 3H sigma_2 - dsigma_1 - dsigma_2
       - sigma_1^2 - sigma_1 sigma_2 - sigma_2^2.
```

Each component is a nonzero polynomial generically.  All components vanish in
the isotropic FRW limit:

```text
sigma_1 = sigma_2 = dsigma_1 = dsigma_2 = 0.
```

## Machine-check export

Execution script:

```text
fundamental_action_reconstruction/p1974_s924_strict_bianchi_anisotropic_eom_residual_obstruction_witness.py
```

Generated witness:

```text
fundamental_action_reconstruction/generated/p1974_s924_strict_bianchi_anisotropic_eom_residual_obstruction_witness.json
```

The JSON exports:

1. the Bianchi-I/FRW component formulas,
2. the symbolic anisotropic residual vector,
3. polynomial nonzero flags,
4. exact isotropic-limit zero residuals,
5. three deterministic rational anisotropic probes with nonzero `scipy` norms,
6. a false-pass guard blocking tensor-level background-independence promotion.

## Theorem statement

The scalar finite-part transport proved in `P1973` is insufficient for full
FRW/Bianchi-I metric EOM covariance.  On a generic diagonal Bianchi-I branch,
the anisotropic Einstein residual vector is polynomially nonzero.  Therefore
background-independence remains obstructed unless the strict theory supplies
one of the following:

1. an anisotropic stress/source tensor `Pi_ij^strict` with zero FRW limit and
   Bianchi-I components cancelling the exported residual vector, or
2. a genuinely tensorial transport connection that maps the full EOM operator
   bundle, not just the scalar finite-part multiplier.

## False-pass boundary

This packet does not prove:

1. global background-independence,
2. PO2 sufficiency / `DELTA_BG_Yf` closure,
3. PO3 nonempty global admissible class,
4. BRST closure,
5. Cutkosky/unitarity closure,
6. `QW-2191` selector discharge,
7. ToE closure.

It proves an obstruction that must be handled before any of those promotions.

## Next honest step

Construct or rule out a strict anisotropic stress/source tensor
`Pi_ij^strict` whose FRW limit is zero and whose diagonal Bianchi-I components
cancel the exported residual vector componentwise.  If no strict-side source
exists, record a no-go theorem for this background-independence route.

## Lay explanation

W poprzednim kroku sprawdziliśmy prosty most między dwoma geometriami i most
zadziałał dla jednego wspólnego mnożnika.  Teraz sprawdziliśmy prawdziwszy test:
konkretne równania grawitacji w geometrii anizotropowej.  Pojawiły się dodatkowe
składniki, których prosty mnożnik nie usuwa.  To nie jest porażka — to uczciwe
wskazanie, czego teoria musi jeszcze dostarczyć, żeby mogła rościć sobie prawo
do pełnej niezależności od tła.
