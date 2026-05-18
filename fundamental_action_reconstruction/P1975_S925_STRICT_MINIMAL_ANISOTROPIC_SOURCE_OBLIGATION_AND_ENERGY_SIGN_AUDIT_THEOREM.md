# P1975 / S925 — Strict Minimal Anisotropic Source Obligation and Energy-Sign Audit Theorem Packet

Status: `OPEN_OBSTRUCTION_WITH_TRACE`  
Route: `strict_only`  
Date anchor: `2026-05-18`

## Professor-level decision

`P1974` showed that scalar FRW/Bianchi-I finite-part transport is insufficient
for tensor-level EOM covariance.  The next honest step is not to invent a new
strict source.  The next honest step is to compute the **minimal source tensor
obligation** that would be required to cancel the `P1974` residual and then test
whether its sign is physically admissible.

## Convention

Use the component equation convention:

```text
G_component = T_component.
```

Therefore, if `P1974` exports

```text
R_component := G_component^BI - G_component^FRW,
```

then the minimal cancelling source must satisfy

```text
Delta T_component = R_component.
```

This is only an obligation target.  It is not a derived source until exported
from strict `K_strict -> coefficients -> L_total -> EOM`.

## Minimal source obligation

Let

```text
Q_shear := sigma_1^2 + sigma_1 sigma_2 + sigma_2^2.
```

From `P1974`, the minimal required source components are:

```text
Delta rho = -Q_shear
Delta p_1 =  3H sigma_1 + dsigma_1 - Q_shear
Delta p_2 =  3H sigma_2 + dsigma_2 - Q_shear
Delta p_3 = -3H sigma_1 - 3H sigma_2 - dsigma_1 - dsigma_2 - Q_shear.
```

The spatial part splits into:

```text
mean pressure shift = -Q_shear
traceless pressure shift =
[3H sigma_1 + dsigma_1,
 3H sigma_2 + dsigma_2,
 -3H sigma_1 - 3H sigma_2 - dsigma_1 - dsigma_2]
```

with exactly zero traceless sum.

## Conditional cancellation

If this source is admitted, then:

```text
R_component - Delta T_component = 0
```

componentwise.  It also vanishes in the FRW isotropic limit:

```text
sigma_1 = sigma_2 = dsigma_1 = dsigma_2 = 0.
```

## Energy-sign obstruction

The required energy component is:

```text
Delta rho = -Q_shear.
```

The quadratic form has matrix:

```text
[[1, 1/2],
 [1/2, 1]]
```

with eigenvalues:

```text
1/2, 3/2.
```

Thus `Q_shear > 0` for every nonzero trace-free shear branch, and hence:

```text
Delta rho < 0
```

for every nonzero branch under this convention.

## Machine-check export

Execution script:

```text
fundamental_action_reconstruction/p1975_s925_strict_minimal_anisotropic_source_obligation_and_energy_sign_audit.py
```

Generated witness:

```text
fundamental_action_reconstruction/generated/p1975_s925_strict_minimal_anisotropic_source_obligation_and_energy_sign_audit.json
```

The witness exports:

1. the minimal source components,
2. the post-source residual `[0,0,0,0]` under conditional admission,
3. the zero FRW isotropic limit,
4. the positive-definite shear quadratic form and eigenvalues,
5. rational numeric probes with negative required energy density,
6. a false-pass guard blocking background-independence closure.

## Theorem statement

`P1975` proves a conditional and a negative statement:

1. **Conditional:** if a strict source exactly equal to the exported minimal
   `Delta T` is derived, then it cancels the `P1974` residual componentwise and
   vanishes in the FRW limit.
2. **Negative:** under `G=T`, that minimal source has energy density
   `Delta rho = -Q_shear`, strictly negative for nonzero trace-free shear.

Therefore no positive-energy anisotropic source can be the minimal cancelling
provider in this convention.  A strict theory continuation must either derive
this sign from an allowed gravitational/shear sector, change the tensorial
transport mechanism, or prove a no-go theorem.

## False-pass boundary

This packet does not prove:

1. a derived strict anisotropic source,
2. positive-energy source closure,
3. global background-independence,
4. PO2/PO3 closure,
5. BRST closure,
6. Cutkosky/unitarity closure,
7. `QW-2191` selector discharge,
8. ToE closure.

## Next honest step

Audit the full strict `L_total` sector list for a derived shear-energy term with
exactly this negative `rho_required` signature.  If no such strict-side term is
available, export a no-go theorem: no positive-energy strict anisotropic
provider can cancel the `P1974` Bianchi-I residual in the current route.

## Lay explanation

Po znalezieniu anizotropowego błędu policzyliśmy, jakie dodatkowe źródło musiałoby
istnieć, żeby ten błąd dokładnie usunąć.  Matematycznie takie źródło jest jasne,
ale ma cenę: jego energia wychodzi ujemna dla każdej niezerowej anizotropii.
To oznacza, że nie wolno go po prostu dopisać jako normalnej materii.  Teoria
musi albo wyprowadzić taki znak z geometrii/grawitacji, albo przyznać, że ta
ścieżka niezależności od tła jest zablokowana.
