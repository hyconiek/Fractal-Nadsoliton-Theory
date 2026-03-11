# T2 Sigma-Int To Residual Datum Bridge Theorem Spec

Status: `T2_PACKET_READY_SIGMA_INT_TO_RESIDUAL_DATUM_BRIDGE_THEOREM_SPEC_NO_FALSE_PASS`
As of: `2026-03-11`

## Goal

After `B4`, `B6`, `C37`, and `C38`, the bridge lane is sharply localized:

- `sigma_int_candidate` exists,
- residual `Z2` slot fit exists on the overlay lane,
- the strict core exports a sign-only export-map object into the residual target
  slot (`F311/N422`),
- but strict-core target-slot population (theta supply) remains absent (`N1/C50`)
  and the current strict sigma-int lane does not derive an actual strict-core
  residual orientation datum (`N7`).

`T2` does not claim that the bridge theorem is already proved.

`T2` does something narrower:
- writes a packet-ready theorem spec for a conditional strict-core bridge
  from `sigma_int_candidate` to the residual orientation datum,
- isolates the minimal assumptions under which such a theorem would make sense,
- separates candidate-fit from theorem-level identification.

## Target theorem

### Informal statement

If the residual orientation datum is shown to be internally representable by the
same strict-core `Z2` object class as `sigma_int_candidate`, and if that
identification is compatible with the selector track rather than only with the
overlay lane, then `sigma_int_candidate` may serve as the internal source of
the residual orientation datum.

### Formal target

```text
T2_SigmaInt_To_ResidualDatum_Bridge_Theorem

Assume:
  A1. a strict-core residual orientation datum slot is explicitly specified,
  A2. a strict-core equivalence or export map identifies that slot with
      sigma_int_candidate,
  A3. the identification is compatible with the selector-track residual
      structure and does not rely only on overlay-lane compatibility,
  A4. the identification respects the anti-overclaim boundary.

Then:
  sigma_int_candidate can serve as a strict-core internal representative of the
  residual orientation datum.
```

## Strict-admissible support

1. `B4`
   - canonical candidate `sigma_int_candidate := chi_FR(gamma_pi1)`
2. `B6`
   - factorized residual `Z2` fit exists on a control lane
3. `B7`
   - overlay compatibility with the mode scaffold
4. `B8`
   - selector-track anti-overclaim audit
5. `C36`
   - axiom-branch to strict-track bridge exists only as overlay
6. `C37`
   - candidate internalization present, strict equivalence absent
7. `C38`
   - theorem-spec/export-spec absent
8. `A10`
   - anti-overclaim boundary

## Minimal lemma DAG

### L1. Sigma-int candidate exists

```text
L1:
The strict track already contains a canonical candidate topological datum
sigma_int_candidate.
```

Support:
- `B4`

### L2. Residual fit exists on an overlay lane

```text
L2:
The residual orientation slot admits a candidate-fit with sigma_int_candidate
on the overlay/control lane.
```

Support:
- `B6`
- `B7`
- `C37`

### L3. Overlay fit is not yet strict-core identification

```text
L3:
Overlay compatibility alone is insufficient for a strict-core bridge theorem.
```

Support:
- `C36`
- `C38`

### L4. A strict-core bridge requires an explicit residual datum slot

```text
L4:
A bridge theorem requires an explicit strict-core target slot or target datum,
not only a control-lane fit.
```

Support:
- `C37`
- `C38`

### L5. A strict-core bridge requires an explicit equivalence or export map

```text
L5:
A bridge theorem requires an explicit equivalence or export map from
sigma_int_candidate to the residual datum; current repo state does not yet
provide such a map.
```

Support:
- `C38`

Update (post-`T148`):
- `F311/N422` exports a sign-only export-map object into the residual target
  slot, but it still does not provide a strict-core theta supply nor an actual
  strict-core residual orientation datum population; so the full bridge theorem
  remains conditional.

## Minimal assumption map

### Physical assumptions

1. the residual orientation datum is genuinely `Z2`-valued,
2. the FR/topological object carried by `sigma_int_candidate` is physically
   relevant for the selector problem,
3. the selector residual is representable without introducing a new unrelated
   datum class.

### Technical assumptions

1. a strict-core target slot is exported,
2. a strict-core equivalence or export map is exported,
3. the map is compatible with the anti-overclaim policy,
4. the bridge does not silently import the axiom-augmented lane as strict-core
   discharge.

## Acceptance skeleton

The theorem spec is acceptable only if all of the following stay explicit:

1. `sigma_int_candidate` is only a candidate object before discharge,
2. overlay compatibility is not counted as strict equivalence,
3. any future bridge must specify both:
   - target datum slot,
   - equivalence/export map,
4. the theorem spec remains conditional until these objects exist,
5. no step is described as discharging `QW-2191`.

## What this theorem would establish if discharged

If discharged, `T2` would establish:
- an internal strict-core representative for the residual orientation datum,
- a precise bridge from topological sigma-int data to the selector residual lane.

It would not by itself establish:
- full gauge uniqueness,
- full theorem-level closure of the selector track,
- full ToE closure,
- resolution of unrelated QFT/RG/gravity global blockages.

## Residual blockers after the spec

Even after `T2` is written, the following remain open:

- strict-core target-slot population (theta supply) for the residual datum,
- strict-core internalization beyond the sign-only export-map object (i.e. an
  actual residual orientation datum population rather than only a residual sign convention),
- discharge of `T1`,
- discharge of `T2`,
- `C32_B2` as a separate negative result about raw overlap export.

## Anti-overclaim

`T2` does not claim:
- theorem-level PASS,
- full-closure PASS,
- that candidate-fit already equals strict equivalence,
- that overlay compatibility already equals internal derivation,
- that `QW-2191` is discharged.

## Product of the step

- theorem-lane packet-ready conditional bridge theorem,
- minimal assumption map,
- minimal lemma DAG,
- explicit separation of candidate-fit from strict-core identification.

## Next step

Natural next move:
- decide whether to pursue:
  - discharge of `T1` as a no-internal-source theorem,
  - or explicit construction of the missing strict-core slot/map objects needed
    by `T2`.
