# P1960 Strict QW-2190 SU(3)/SU(2) Structure Constants and Jacobi Certificate Packet

Status: `PARTIAL_STRICT_ALGEBRA_DATA_EXPORT__BRST_STILL_OPEN`
As of: `2026-05-17`

## Purpose

`P1959` correctly refused to promote the `P1958` local Abelian ghost seed to a
full non-Abelian BRST theorem.  However, it checked only the FAR-local path for
the `QW-2190` source report:

```text
fundamental_action_reconstruction/report_qw2190_kernel_mode_representation_emergence_gate.json
```

The actual source used by `P452/P454` is present at repository root:

```text
report_qw2190_kernel_mode_representation_emergence_gate.json
```

`P1960` therefore performs the next narrow correction: ingest the repo-root
`QW-2190` report and export exact structure constants plus exact Jacobi
certificates for the declared `SU(3) x SU(2) x U(1)` direct-product algebra.

This is an algebra-data closure only.  It is not a full BRST/BV theorem.

## Pre-Execution Grep

Before execution the FAR tree was searched in English and Polish for:

```text
P1960, StructureConstants_SU3_SU2_U1_strict_v1,
JacobiCertificate_SU3_SU2_U1_strict_v1, structure constants,
Jacobi identity, QW-2190 source report, non-Abelian BRST,
stale struktury, tozsamosc Jacobiego, nieabelowy BRST
```

Existing relevant sources:

```text
P1959: obstruction audit for non-Abelian BRST data.
N493/N494/N495: QW-2191 rotations/signs are gauge/conjugation only for QW-2190 embedding audits.
P452/P454: executable sign/O(2) gauge-equivalence audits.
report_qw2190_kernel_mode_representation_emergence_gate.json: present at repo root.
```

`TOE_FINAL_DOCUMENTATION_RELEASE_8_STRICT_FULL.pdf` was ignored.

## Algebra Export

The executor uses the Hermitian convention:

```text
T_a = lambda_a/2      for SU(3)
t_i = sigma_i/2      for SU(2)
2 Tr(T_a T_b) = delta_ab
[T_a,T_b] = i f_abc T_c
```

It exports:

```text
StructureConstants_SU3_SU2_U1_strict_v1
JacobiCertificate_SU3_SU2_U1_strict_v1
```

with:

```text
SU(3): a,b,c in 1..8
SU(2): i,j,k in 1..3
U(1): all structure constants zero
cross-factor constants: zero by direct product
```

## Machine Checks

The exact symbolic checks give:

```text
SU(3) Jacobi checked components: 4096
SU(3) nonzero Jacobi components: 0

SU(2) Jacobi checked components: 81
SU(2) nonzero Jacobi components: 0

exact_algebra_pass = true
```

The QW-2190 embedding numerical audit gives:

```text
b3_orthonormal_residual = 2.709517654683241e-16
b2_orthonormal_residual = 3.6042185851109777e-16
kernel_invariance_residual_su3 = 6.003302281971643e-16
kernel_invariance_residual_su2 = 3.831284866810958e-16
embedded_su3_closure_residual = 2.0077411267522386e-16
embedded_su2_closure_residual = 1.7375088326559935e-16
embedded_su3_su2_cross_commutator_residual = 1.184514349962965e-16
tolerance = 1e-12
```

Therefore:

```text
QW2190_ROOT_REPORT = true
F_TABLE = true
JACOBI = true
EMBEDDING = true
```

## BRST Readiness After P1960

The post-P1960 readiness formula is:

```text
F_TABLE
& JACOBI
& QW2190_REPORT
& EMBEDDING
& CONNECTION_RULES
& BRST_RULES
& GHOST_SELF
& BV_MAP
& GHOST_CUT
```

The truth assignment is:

```text
F_TABLE = true
JACOBI = true
QW2190_REPORT = true
EMBEDDING = true
CONNECTION_RULES = false
BRST_RULES = false
GHOST_SELF = false
BV_MAP = false
GHOST_CUT = false
```

So:

```text
evaluated_brst_ready = false
```

## Scope Guard

`P1960` may claim:

```text
repo-root QW-2190 source report is present and ingested
exact SU(3), SU(2), U(1) direct-product structure constants are exported
exact Jacobi certificates for those constants are exported
QW-2190 embedding numerical residuals pass within declared tolerance
```

`P1960` must not claim:

```text
global QW-2191 discharge
strict-core selector closure
full SU(3)xSU(2)xU(1) BV/BRST operator map
Q^2=0 on the complete strict field bundle
ghost-cancelled Cutkosky theorem
TG2_BRST_GLOBAL_NILPOTENCY PASS
TG3_CUTKOSKY_GLOBAL_UNITARITY PASS
ToE closure
```

## Outputs

- `p1960_s910_strict_qw2190_su3_su2_structure_constants_and_jacobi_certificate.py`
- `generated/p1960_s910_strict_qw2190_su3_su2_structure_constants_and_jacobi_certificate.json`

## Next Honest Step

Build `P1961` with high reasoning: construct the field-level non-Abelian BRST
differential

```text
s A_mu^a, s c^a, s cbar^a, s B^a
```

from the `P1960` constants and the `P1958` gauge-fixing seed, then execute an
explicit graded nilpotency check.  If the strict `L_total` field registry cannot
support that map, export the precise field-registry obstruction.

## Lay Explanation

Znaleziono brakujacy raport `QW-2190` w katalogu glownym repo i zbudowano
dokladna tabele regul komutatorow `SU(3)`, `SU(2)` oraz zerowa czesc `U(1)`.
Sprawdzono tez tozsamosc Jacobiego, czyli warunek, bez ktorego pelne duchy
BRST nie moga dzialac.  To nadal nie jest pelny BRST: mamy teraz algebre, ale
nie cala maszyne dzialania na polach.
