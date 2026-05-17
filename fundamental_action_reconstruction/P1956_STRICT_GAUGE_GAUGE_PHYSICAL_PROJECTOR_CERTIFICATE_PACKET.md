# P1956 Strict Gauge-Gauge Physical Projector Certificate Packet

Status: `OPEN_OBSTRUCTION_WITH_TRACE__TRANSVERSE_PROJECTOR_EXPORTED_NO_BRST_PASS`
As of: `2026-05-17`

## Pre-Execution Grep

Before execution, the repository was searched in English and Polish for:

```text
BRST, physical-state projector, polarization projector, transverse polarization,
ghost sector, ghost cancellation, gauge_gauge, Cutkosky,
projektor, polaryzacja, poprzeczny, sektor duchow, kasowanie duchow,
unitarnosc, kohomologia
```

The search found gate and template sources, but no exported
`BRSTPhysicalProjector_gauge_gauge_strict_B1_v1` and no channel-level
ghost-subtracted BRST cohomology projector:

```text
P1767: BRST/Cutkosky requires G_BW PASS, ghost_sector_nonproxy_export,
       and BV_BRST_operator_map.
P1801: BRST nilpotency intake template, not a channel projector.
P1802: Cutkosky unitarity intake template, not a computed projector.
P1854: B1 BRST cochain and first Cutkosky channel proxy, still pre-theorem.
P1954: M2 BRST projector listed as missing.
P1955: minimal tree-level hAA vertex exported; BRST projection still open.
```

## Result

`P1956` exports the exact local on-shell transverse polarization projector for
the two gauge bosons in the `graviton -> gauge_gauge` cut.

The flat cut frame is:

```text
eta_mu_nu = diag(-1,1,1,1)
k1 = (1,0,0, 1),  n1 = (1,0,0,-1)
k2 = (1,0,0,-1),  n2 = (1,0,0, 1)
eps_x = (0,1,0,0), eps_y = (0,0,1,0)
```

For each gauge boson, the exported projector is:

```text
P^{mu nu}(k,n)
  = eta^{mu nu}
    - (k^mu n^nu + n^mu k^nu)/(k.n)
    + (n.n) k^mu k^nu/(k.n)^2.
```

In the selected lightlike reference frame `n.n=0`, `k.n=-2`, and the exact
matrix becomes:

```text
P^{mu nu} = diag(0,1,1,0).
```

The symbolic check proves:

```text
sum_{lambda=x,y} eps_lambda^mu eps_lambda^nu - P^{mu nu}(k,n) = 0
P^{mu nu} k_nu = 0
P^{mu nu} n_nu = 0
P^mu_alpha P^alpha_nu - P^mu_nu = 0
rank(P)=2
```

For the two-particle final state:

```text
P_gauge_gauge = P_1 mixed kron P_2 mixed
rank(P_gauge_gauge)=4
trace(P_gauge_gauge)=4
P_gauge_gauge^2 - P_gauge_gauge = 0
```

## Scope

This discharges the `P1954` missing `M2` item only in the limited sense:

```text
M2_BRST_projector =
  DISCHARGED_ONLY_FOR_LOCAL_ONSHELL_TRANSVERSE_POLARIZATION_SUM
```

It does not provide:

1. a BV/BRST charge `Q`,
2. a `Q^2=0` nilpotency proof,
3. a nonproxy ghost-sector action,
4. a ghost-cancellation trace,
5. a proof that the full cut sum is restricted to BRST cohomology,
6. a dressed `DiscM=CutSum` theorem,
7. a global `UR_link`.

Therefore `TG2_BRST_GLOBAL_NILPOTENCY` and
`TG3_CUTKOSKY_GLOBAL_UNITARITY` are explicitly **not promoted**.

## Outputs

- `p1956_s906_strict_gauge_gauge_physical_projector_certificate.py`
- `generated/p1956_s906_strict_gauge_gauge_physical_projector_certificate.json`

## Next Honest Step

Build `P1957` with high reasoning: either construct the strict BV/BRST charge,
ghost-sector nonproxy action, `Q^2=0` nilpotency check, and
`gauge_gauge` ghost-cancellation trace, or export a formal obstruction proving
that those data are not yet available.
