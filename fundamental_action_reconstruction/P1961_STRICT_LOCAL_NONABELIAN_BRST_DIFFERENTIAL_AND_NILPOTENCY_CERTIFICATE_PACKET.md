# P1961 Strict Local Non-Abelian BRST Differential and Nilpotency Certificate Packet

Status: `LOCAL_NONABELIAN_BRST_DIFFERENTIAL_EXPORTED__GLOBAL_BV_BRST_STILL_OPEN`
As of: `2026-05-17`

## Purpose

`P1960` exported exact `SU(3)`, `SU(2)`, and `U(1)` structure constants plus
Jacobi certificates.  `P1961` uses those data to construct the local
gauge-sector BRST differential:

```text
s A_mu^a = partial_mu c^a + g f^a_bc A_mu^b c^c
s c^a    = -g/2 f^a_bc c^b c^c
s cbar^a = B^a
s B^a    = 0
```

and checks `s^2 = 0` on the local generators:

```text
A_mu^a, c^a, cbar^a, B^a
```

for the direct-product `SU(3) x SU(2) x U(1)` gauge sector.

This is a local gauge-sector theorem only.  It is not a full BV theorem, not a
global BRST charge construction, not a cohomology theorem, and not a Cutkosky
unitarity closure.

## Pre-Execution Grep

Before execution, the FAR tree was searched in English and Polish for:

```text
P1961, NonAbelianBRSTDifferential_strict_B1_v1,
GaugeConnectionRules_SU3_SU2_U1_strict_v1,
GhostSelfInteraction_strict_B1_v1, LocalNilpotencyCertificate,
graded Leibniz, BRST differential, nilpotency,
rozniczka BRST, nilpotencja, akcja duchow
```

Relevant existing sources:

```text
P1958: local Abelian Lorenz/FP ghost seed.
P1959: non-Abelian BRST was blocked before structure constants/Jacobi.
P1960: exact SU(3), SU(2), U(1) structure constants and Jacobi certificates.
```

`TOE_FINAL_DOCUMENTATION_RELEASE_8_STRICT_FULL.pdf` was ignored.

## Exported Local Rules

`P1961` exports:

```text
GaugeConnectionRules_SU3_SU2_U1_strict_v1
NonAbelianBRSTDifferential_strict_B1_v1
GhostSelfInteraction_strict_B1_v1
LocalNilpotencyCertificate_Q2_zero_gauge_sector_strict_B1_v1
```

The local Faddeev-Popov operator is recorded as:

```text
M_FP^a_c = partial^mu( delta^a_c partial_mu + g f^a_bc A_mu^b )
```

and the local ghost Lagrangian is:

```text
L_ghost = - cbar_a partial^mu( partial_mu c^a + g f^a_bc A_mu^b c^c )
```

For `U(1)`, `f = 0`, so this reduces to the `P1958` Abelian seed:

```text
L_ghost = - cbar Box c
```

## Machine Certificate

The executor uses a coefficient-level exterior-sign engine.  It canonicalizes
odd ghost monomials and checks the BRST square coefficient by coefficient.

Results:

```text
SU(3): all local generator s^2 checks pass
SU(2): all local generator s^2 checks pass
U(1): all local generator s^2 checks pass

all_local_generator_s2_zero = true
evaluated_local_gauge_sector_ready = true
evaluated_global_TG2_ready = false
```

For `SU(3)`, the largest checked blocks were:

```text
s2_c_SU3_all: checked_targets = 8
s2_A_derivative_terms_SU3_all_mu: checked_targets = 32
s2_A_connection_terms_SU3_all_mu: checked_targets = 256
```

All nonzero residual samples are empty.

## Scope Guard

`P1961` may claim:

```text
local non-Abelian gauge-sector BRST differential exported
formal local FP ghost self-interaction exported
exact s^2=0 on local gauge-sector generators A,c,cbar,B
U(1) reduction matches the P1958 Abelian seed
```

`P1961` must not claim:

```text
global BV/BRST theorem
BRST charge Q on the full strict phase space
cohomology or physical-state projection
full L_total gauge invariance
ghost-cancelled Cutkosky theorem
TG2_BRST_GLOBAL_NILPOTENCY PASS
TG3_CUTKOSKY_GLOBAL_UNITARITY PASS
QW-2191 global selector discharge
ToE closure
```

## Outputs

- `p1961_s911_strict_local_nonabelian_brst_differential_and_nilpotency_certificate.py`
- `generated/p1961_s911_strict_local_nonabelian_brst_differential_and_nilpotency_certificate.json`

## Next Honest Step

Build `P1962` with high reasoning: audit the strict `L_total` field registry
and matter/Higgs/spinor representations to see whether the `P1961` local
gauge-sector BRST differential extends to the full exported nonproxy bundle. If
it cannot, export the precise representation/field-registry obstruction.

## Lay Explanation

Po `P1960` mielismy tablice regul symetrii.  `P1961` pokazuje, ze lokalna
maszyna BRST dla samych pol cechowania dziala algebraicznie: dwa kolejne kroki
BRST daja zero.  To jeszcze nie jest pelny dowod dla calej teorii, bo trzeba
sprawdzic, jak te reguly dzialaja na materie, Higgs/spinory, grawitacje i pelny
`L_total`.
