# T10 Route-Role Typing Rule Theorem Spec

Status: `T10_PACKET_READY_ROUTE_ROLE_TYPING_RULE_THEOREM_SPEC_NO_FALSE_PASS`
As of: `2026-03-06`

## Goal

`T9` reduced the failure of `T8` to one meta-level blocker:

- `T9_B1 := no formal route-role typing rule or admissibility-by-role declaration showing that every current strict-core theta-export route must instantiate exactly one of the six named route roles`

`T10` does not claim that such a typing rule already exists.

`T10` does something narrower:
- writes a packet-ready theorem spec for the missing route-role typing rule,
- isolates the minimal typing clauses required to turn the present role vocabulary into an admissibility-by-role declaration,
- makes explicit how such a rule would discharge `T8` and therefore feed back into `T6`, `T4`, and `T1`.

## Target theorem

### Informal statement

For the current strict-core selector track, every admissible current theta-export
route has exactly one route-role type among:

1. raw overlap,
2. formula,
3. representative,
4. downstream schema,
5. source skeleton,
6. strict-to-axiom bridge.

Hence admissibility is exhausted by those six role types.

### Formal target

```text
T10_RouteRole_Typing_Rule_Theorem

Assume:
  A1. the present selector track exposes a finite role vocabulary
      {raw_overlap, formula, representative, downstream_schema,
       source_skeleton, strict_to_axiom_bridge},
  A2. every current admissible theta-export route carries exactly one role type
      from that vocabulary,
  A3. no additional role type is exported elsewhere in the present strict core,
  A4. admissibility-by-role is evaluated under the current anti-overclaim boundary.

Then:
  admissibility of current strict-core theta-export routes is exhausted by the
  six named route-role types.
```

## Strict-admissible support

1. `C32`
   - raw overlap route explicitly isolated
2. `C33`
   - formula-class route explicitly isolated
3. `C34`
   - representative-class route explicitly isolated
4. `C49`
   - downstream populated-instance route explicitly isolated
5. `C50`
   - strict-core source-skeleton route explicitly isolated as absent
6. `C51`
   - strict-to-axiom bridge route explicitly isolated as fallback-only
7. `T8`
   - admissibility grammar requires a route-role basis
8. `T9`
   - discharge attempt isolates the missing route-role typing rule
9. `A10`
   - anti-overclaim boundary

## Minimal lemma DAG

### L1. Present route-role vocabulary is finite and explicit

```text
L1:
The current selector track already exports a finite explicit route-role
vocabulary consisting of exactly six named route-role types.
```

Support:
- `C32`
- `C33`
- `C34`
- `C49`
- `C50`
- `C51`

### L2. Each audited route instance carries one stable role label

```text
L2:
Each currently audited theta-export route instance is stably classified by one
role label from the six-role vocabulary.
```

Support:
- `C32`
- `C33`
- `C34`
- `C49`
- `C50`
- `C51`

### L3. No additional role type is currently exported in strict core

```text
L3:
No extra route-role type for actual theta-export is currently exported outside
those six named role labels.
```

Support:
- `T9`
- `A10`

### L4. Role-typing closure upgrades vocabulary to admissibility-by-role

```text
L4:
If L1-L3 hold, then every current admissible theta-export route is typed by the
six-role vocabulary, so admissibility is exhausted by role type.
```

Support:
- `T8`
- `T9`

## Minimal assumption map

### Technical assumptions

1. route-role labels are intrinsic to route construction, not post hoc tags,
2. each admissible route has exactly one role type in the current selector track,
3. no hidden role type is introduced by reinterpretation of downstream or fallback lanes,
4. admissibility ranges only over routes expressible in the present strict-core selector formalism.

### Anti-overclaim assumptions

1. the theorem ranges only over the present repo state,
2. the theorem does not quantify over future extensions,
3. the theorem does not convert the axiom-augmented lane into strict-core discharge,
4. the theorem does not discharge `QW-2191` by itself.

## Acceptance skeleton

The theorem spec is acceptable only if all of the following stay explicit:

1. route-role typing is defined only for the present selector track,
2. the six role labels are treated as current candidates, not eternal universals,
3. hidden-role exclusion is restricted to current exported strict-core objects,
4. fallback citation remains non-strict,
5. no theorem-level or full-closure PASS language is introduced.

## What this theorem would establish if discharged

If discharged, `T10` would establish:
- a formal route-role typing rule for current theta-export routes,
- the missing admissibility-by-role step required by `T9`,
- a credible route to discharging `T8`, then `T6`, then `T4`, then `T1`.

It would not by itself establish:
- `T2`,
- full selector-track closure,
- full gauge uniqueness,
- full ToE closure.

## Residual blockers after the spec

Even after `T10` is written, the following remain open:
- actual discharge of `T10`,
- actual discharge of `T8`,
- actual discharge of `T6`,
- actual discharge of `T4`,
- actual discharge of `T1`,
- actual discharge of `T2`,
- `C32_B2` as an independent negative result about the raw overlap route.

## Anti-overclaim

`T10` does not claim:
- theorem-level PASS,
- full-closure PASS,
- that the route-role typing rule is already proved,
- that `T8`, `T6`, `T4`, or `T1` are already discharged,
- that `QW-2191` is discharged.

## Product of the step

- packet-ready theorem spec for the missing route-role typing rule,
- minimal lemma DAG for admissibility-by-role closure,
- explicit acceptance skeleton,
- maintained no-false-pass discipline.

## Next step

Natural next move:
- perform a first discharge attempt for `T10` by checking whether the current
  selector-track audits already imply uniqueness and exhaustiveness of route-role typing.
