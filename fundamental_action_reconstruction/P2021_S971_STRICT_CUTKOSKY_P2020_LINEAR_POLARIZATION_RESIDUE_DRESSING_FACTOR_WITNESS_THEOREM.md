# P2021 S971 Strict Cutkosky P2020 Proxy Residue Admissibility Audit

Status: `P2021_PROXY_RESIDUE_TRANSPORT_SANITY_ONLY_NOT_P1953_ADMISSIBLE_NO_FALSE_PASS`
As of: `2026-05-18`

## Goal

P2020 exported the exact tree-level, phase-space-integrated, real
linear-polarization-resolved CutSum matrix for

```text
gravity/gauge channel: graviton -> gauge_gauge.
```

The reviewed P2021 correction is deliberately stricter than the first draft.
P1994's

```text
Z(s) = (51*s + 50)/(50*(s + 1))
```

is a proxy/mirrored import, not a loop-derived same-scheme dressed residue.
P2021 therefore keeps the scalar-transport algebra only as a sandbox and
**rejects** the proxy as a valid `P1953` dressed-interface update.

## P2020 matrix input

In units of `kappa^2*Z_gauge^2`, P2020 supplies

```text
C_no_sym = [[1/pi,0],[0,1/pi]],
C_ident  = [[1/(2*pi),0],[0,1/(2*pi)]],
```

where the basis order is the real graviton linear-polarization basis

```text
{plus, cross}.
```

The scalar values `2/pi` and `1/pi` remain traces of these matrices, not the
primary theorem object.

## Local algebra sandbox

For real `s > 0`, the proxy factor satisfies

```text
Z_proxy(s) = (51*s + 50)/(50*(s + 1)) > 0,
R_proxy(s) = Z_proxy(s)^2 > 0.
```

Multiplying the P2020 positive matrix by `R_proxy(s)` preserves positivity and
does not rotate the linear-polarization basis:

```text
C_no_sym_proxy(s)
  = R_proxy(s) C_no_sym
  = [[(51*s + 50)^2/(2500*pi*(s + 1)^2), 0],
     [0, (51*s + 50)^2/(2500*pi*(s + 1)^2)]],

C_ident_proxy(s)
  = R_proxy(s) C_ident
  = [[(51*s + 50)^2/(5000*pi*(s + 1)^2), 0],
     [0, (51*s + 50)^2/(5000*pi*(s + 1)^2)]].
```

This statement is intentionally weak: it is a matrix-preserving algebra check,
not a dressed Cutkosky theorem.

## Admissibility failure

P2021 exports the following negative admissibility result for the P1953 dressed
interface:

```text
loop_derived_from_L_total = false
same_scheme_as_P1953_MSbar_B1_seed_locked = false
DiscM_common_basis_exported = false
DiscM_minus_CutSum_simplified_evaluated = false
BRST_physical_state_projector_exported = false
proxy_not_promoted_to_dressed_residue = true
```

Therefore the proxy transport is **not** accepted as any of the missing P1953
objects:

```text
M_dressed_common_basis
  = OPEN_NOT_EXPORTED__P2021_PROXY_FACTOR_REJECTED_AS_DRESSED_INPUT
AbsM_dressed_squared_common_basis
  = OPEN_NOT_EXPORTED__PROXY_TRANSPORT_SANDBOX_ONLY
CutSum_common_basis
  = P2020_TREE_LINEAR_POLARIZATION_CUTSUM_REMAINS_AVAILABLE_BUT_NOT_DRESSED
DiscM_common_basis
  = OPEN_NOT_EVALUATED
DiscM_minus_CutSum_simplified
  = OPEN_NOT_EVALUATED
```

## Numerical certificate

The script evaluates the proxy-transported no-identical-symmetry matrix on the
bounded sample grid

```text
s in {1/2, 1, 2, 4, 8}
```

and checks that the numerical eigenvalues agree with the exact symbolic
eigenvalues to below `1e-14` in L2 norm.  This verifies only the sandbox algebra.
It does not verify loop unitarity.

## No-false-pass boundary

P2021 remains an obstruction-trace audit.  It proves only that a positive scalar
proxy factor would preserve the P2020 matrix if such a factor were admissible.
It simultaneously proves that the current P1994 proxy is **not** admissible as a
P1953 dressed residue.

P2021 is not:

1. a loop-derived dressed amplitude,
2. a `DiscM_common_basis` calculation,
3. a BRST cohomology projector,
4. a same-scheme `DiscM = CutSum` theorem,
5. all-state unitarity,
6. `QW-2191` discharge,
7. ToE closure.

## Progress toward ToE

This is a correction in rigor rather than a larger physics result.  It prevents
a convenient positive proxy factor from being mistaken for the missing dressed
unitarity data.  The useful local algebra survives, but the theorem frontier is
now sharper: the next object must be loop-derived and scheme-locked, not a
mirrored scalar proxy.

## Next honest step

Build P2022 as a real same-scheme source audit/derivation: either extract a
loop-derived `hAA` self-energy/discontinuity residue from `L_total` with BRST
projector data, or export a nonavailability theorem listing the exact missing
vertices, gauges, and projector inputs needed before `DiscM_minus_CutSum` can be
evaluated.
