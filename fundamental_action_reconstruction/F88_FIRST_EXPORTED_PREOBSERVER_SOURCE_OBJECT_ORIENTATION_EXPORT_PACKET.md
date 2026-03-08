# F88 First Exported Preobserver Source Object Orientation Export Packet

Status: `F88_EXECUTED_FIRST_EXPORTED_PREOBSERVER_SOURCE_OBJECT_ORIENTATION_EXPORT_PACKET_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

After `N195`, the source-object side is no longer the main blocker on this
lane.

The next honest constructive move is:

```text
export one explicit preobserver orientation datum from
S_preLM_strict_core_source_object_v1
without pretending that B_sel, R_sel, O_sel, or selector closure already exist
```

## Source object reused

Reuse:

```text
S_preLM_strict_core_source_object_v1
  = u_T + cos(phi) u_L + (cos(phi)/4) u_M
```

on the typed carrier:

```text
V_topo ⊕ L_int ⊕ M_int
```

with future orientation-export target frame from `F82`:

```text
E_orient_target_frame_v1 := (u_T, u_L)
```

## Topological-light projection

Project onto the strict preobserver topological-light plane:

```text
pi_TL(S_preLM_strict_core_source_object_v1) = u_T + cos(phi) u_L
```

Let:

```text
N_TL := sqrt(1 + cos(phi)^2)
```

## Exported orientation datum

Freeze one explicit preobserver orientation export:

```text
e_parallel := (u_T + cos(phi) u_L) / N_TL
e_transverse := (-cos(phi) u_T + u_L) / N_TL
```

and define:

```text
E_orient_preLM_v1 := (e_parallel, e_transverse)
```

This is an ordered orthonormal frame on `span{u_T, u_L}` derived from the
exported source object.

## Why this is an honest internal export

1. it is derived only from `S_preLM_strict_core_source_object_v1`,
2. it uses no observer slot,
3. it uses no imported `psi0`,
4. it uses no external selector control,
5. it uses only the strict preobserver carrier declared in `F82`,
6. it keeps `K_strict_gate` at operational-coefficient scope only.

## Quotient / gauge-safety witness

Under any positive rescaling `lambda > 0` of the topological-light projection:

```text
normalize(lambda * (u_T + cos(phi) u_L)) = e_parallel
```

So the exported direction is invariant under positive scalar rescaling of the
upstream topological-light support.

## Bridge-ready output

Freeze the immediate bridge-ready structure:

```text
B_sel_start_frame_v1 := E_orient_preLM_v1
selector_axis_v1 := e_parallel
```

This still does not export `B_sel`.

It only makes the next bridge start typed and explicit.

## Hard limits

`F88` does not claim:

- actual `B_sel`,
- actual `R_sel`,
- actual `O_sel`,
- downstream completion,
- strict-core selector closure,
- `QW-2191` discharge,
- ToE closure.

## Recommended next move

The correct next move is:

```text
test the full F32 admissibility contract directly on E_orient_preLM_v1
```
