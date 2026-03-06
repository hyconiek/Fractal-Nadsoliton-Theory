# T12 Typing Judgment Totality and Uniqueness Theorem Spec

Status: `T12_PACKET_READY_TYPING_JUDGMENT_TOTALITY_UNIQUENESS_THEOREM_SPEC_NO_FALSE_PASS`
As of: `2026-03-06`

## Goal

`T11` reduced the failure of `T10` to one meta-level blocker:

- `T11_B1 := no formal typing judgment or totality-and-uniqueness clause showing that every current admissible strict-core theta-export route has exactly one route-role label in the six-role vocabulary`

`T12` does not claim that such a typing judgment already exists.

`T12` does something narrower:
- writes a packet-ready theorem spec for the missing formal typing judgment,
- isolates the minimal totality and uniqueness clauses required to turn the present role vocabulary into a theorem-level route-role typing rule,
- makes explicit how such a judgment would discharge `T10` and therefore feed back into `T8`, `T6`, `T4`, and `T1`.

## Target theorem

### Informal statement

For the current strict-core selector track there exists a formal judgment

```text
route ⊢ role
```

such that:
- every current admissible theta-export route has at least one role label from
  the six-role vocabulary,
- no current admissible theta-export route has more than one such role label.

Hence every current admissible route has exactly one role label in the present
six-role vocabulary.

### Formal target

```text
T12_Typing_Judgment_Totality_Uniqueness_Theorem

Assume:
  A1. there exists a formal judgment form `route ⊢ role` for the present
      selector track,
  A2. (totality) every current admissible theta-export route satisfies
      `route ⊢ role` for at least one role in
      {raw_overlap, formula, representative, downstream_schema,
       source_skeleton, strict_to_axiom_bridge},
  A3. (uniqueness) no current admissible theta-export route satisfies
      `route ⊢ role_1` and `route ⊢ role_2` with `role_1 != role_2`,
  A4. no extra role outside the six-role vocabulary is exported in the current
      strict core,
  A5. admissibility is evaluated under the current anti-overclaim boundary.

Then:
  every current admissible strict-core theta-export route has exactly one
  route-role label in the six-role vocabulary.
```

## Strict-admissible support

1. `C32`
   - raw overlap route explicitly isolated and stably labelled
2. `C33`
   - formula-class route explicitly isolated and stably labelled
3. `C34`
   - representative-class route explicitly isolated and stably labelled
4. `C49`
   - downstream populated-instance route explicitly isolated and stably labelled
5. `C50`
   - source-skeleton route explicitly isolated as absent
6. `C51`
   - strict-to-axiom bridge route explicitly isolated as fallback-only
7. `T10`
   - route-role typing rule requires exact role typing
8. `T11`
   - discharge attempt isolates the missing formal typing judgment with totality and uniqueness
9. `A10`
   - anti-overclaim boundary

## Minimal lemma DAG

### L1. Present role vocabulary is finite and explicit

```text
L1:
The current selector track already exports a finite explicit six-role vocabulary.
```

Support:
- `C32`
- `C33`
- `C34`
- `C49`
- `C50`
- `C51`

### L2. Known route instances exhibit stable role labels

```text
L2:
Each currently audited route instance is stably associated with one role label
from the six-role vocabulary.
```

Support:
- `C32`
- `C33`
- `C34`
- `C49`
- `C50`
- `C51`

### L3. Hidden extra role labels are excluded in current strict core

```text
L3:
No extra role label outside the six-role vocabulary is currently exported.
```

Support:
- `T11`
- `A10`

### L4. Typing judgment + totality + uniqueness upgrades labels to exact role typing

```text
L4:
If L1-L3 hold and a formal judgment `route ⊢ role` satisfies totality and
uniqueness, then every current admissible theta-export route has exactly one
role label in the six-role vocabulary.
```

Support:
- `T10`
- `T11`

## Minimal assumption map

### Technical assumptions

1. a formal judgment form `route ⊢ role` exists for the present selector track,
2. totality ranges over all current admissible theta-export routes,
3. uniqueness excludes multi-typing inside the six-role vocabulary,
4. no hidden extra role is introduced by reinterpretation of downstream or fallback lanes.

### Anti-overclaim assumptions

1. the theorem ranges only over the present repo state,
2. the theorem does not quantify over future extensions,
3. the theorem does not convert the axiom-augmented lane into strict-core discharge,
4. the theorem does not discharge `QW-2191` by itself.

## Acceptance skeleton

The theorem spec is acceptable only if all of the following stay explicit:

1. the typing judgment ranges only over the present selector track,
2. totality and uniqueness are stated explicitly, not implied rhetorically,
3. hidden-role exclusion is restricted to current exported strict-core objects,
4. fallback citation remains non-strict,
5. no theorem-level or full-closure PASS language is introduced.

## What this theorem would establish if discharged

If discharged, `T12` would establish:
- a formal route-role typing judgment for current theta-export routes,
- explicit totality and uniqueness clauses,
- the missing exact-typing step required by `T11`,
- a credible route to discharging `T10`, then `T8`, then `T6`, then `T4`, then `T1`.

It would not by itself establish:
- `T2`,
- full selector-track closure,
- full gauge uniqueness,
- full ToE closure.

## Residual blockers after the spec

Even after `T12` is written, the following remain open:
- actual discharge of `T12`,
- actual discharge of `T10`,
- actual discharge of `T8`,
- actual discharge of `T6`,
- actual discharge of `T4`,
- actual discharge of `T1`,
- actual discharge of `T2`,
- `C32_B2` as an independent negative result about the raw overlap route.

## Anti-overclaim

`T12` does not claim:
- theorem-level PASS,
- full-closure PASS,
- that the typing judgment is already proved,
- that `T10`, `T8`, `T6`, `T4`, or `T1` are already discharged,
- that `QW-2191` is discharged.

## Product of the step

- packet-ready theorem spec for the missing formal typing judgment,
- explicit totality and uniqueness clauses,
- minimal lemma DAG for exact route-role typing,
- maintained no-false-pass discipline.

## Next step

Natural next move:
- perform a first discharge attempt for `T12` by checking whether the current
  route audits already imply a usable judgment form plus totality and uniqueness.
