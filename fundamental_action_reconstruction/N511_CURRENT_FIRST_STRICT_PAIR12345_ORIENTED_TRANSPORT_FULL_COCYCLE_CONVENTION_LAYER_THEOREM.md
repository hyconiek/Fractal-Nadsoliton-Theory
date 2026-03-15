# N511 Current First Strict `pair1..pair5` Oriented Transport (α mod 2π) — Full Cocycle Convention-Layer Theorem (No False‑PASS)

Status: `N511_DISCHARGED_CURRENT_FIRST_STRICT_PAIR12345_ORIENTED_TRANSPORT_FULL_COCYCLE_CONVENTION_LAYER_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`F466` exports a lane‑scoped five‑chart selector‑atlas ingredient on `{pair1..pair5}` at projector level with full triple cocycle/path‑independence data
(packaged by `N510`). This is axis-only (angles `mod π`) and therefore sign‑gauge‑safe, but it does not supply oriented transport (`mod 2π`) at the level of
representative vectors.

`F467` performs the next honest strict move:

```text
export a tracked gauge/convention lift of the five‑chart atlas transport to oriented angles α mod 2π,
induced by the currently exported representative vectors u_1..u_5,
and record full triple cocycle/path‑independence at vector level.
```

This theorem packages that statement in strict no‑false‑PASS discipline.

## Strict‑admissible inputs reused

1. `F467`
   - exported oriented theta family (mod 2π; tracked convention),
   - exported oriented transport operators on `{pair1..pair5}`,
   - exported oriented vector section with full triple cocycle audits.
2. `P470`
   - probe‑level audit that the exported oriented transport operators satisfy vector transport and full triple cocycle/path‑independence relations
     (numeric tolerance on the current exported instance).
3. `N502`
   - residual sign can be frozen as a tracked gauge/convention layer where sign is gauge-irrelevant; this theorem is used only as discipline:
     **the oriented lift remains a convention layer and is not promoted to a physical sign-sensitive datum**.
4. `A10`
   - anti‑overclaim boundary.

## Theorem (oriented transport lift exists as a tracked convention layer)

From `F467`, the repo exports:

1. an explicit oriented theta family lift `θ_m mod 2π` on `{pair1..pair5}` induced by the currently exported representative vectors `u_1..u_5`,
2. an explicit family of oriented chart transport operators (`α mod 2π`) on `{pair1..pair5}`,
3. a vector section packaging the exported `u_1..u_5` together with vector-level transport audits, and
4. explicit full triple cocycle/path‑independence audit data at the level of those exported vectors.

Probe-level evidence that the exported current instance satisfies the declared vector transport and cocycle equations is recorded by `P470`.

Therefore the strict core now contains a lane‑scoped oriented-transport lift on `{pair1..pair5}` as a **tracked gauge/convention layer** with full triple cocycle data at vector level. ∎

## What `N511` does not claim

`N511` does not claim:

1. a sign‑sensitive physical orientation datum (lifting residual `Z2` as physics),
2. a global selector atlas on the full strict domain `C_v1`,
3. strict-core selector closure / admissible `S_sel_int`,
4. global discharge of `QW-2191`,
5. ToE closure.

