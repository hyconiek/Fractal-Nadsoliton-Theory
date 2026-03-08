# T14 Source Topology Selector Theorem Spec

Status: `T14_PACKET_READY_SOURCE_TOPOLOGY_SELECTOR_THEOREM_SPEC_FUTURE_ROUTE_ONLY_NO_FALSE_PASS`
As of: `2026-03-08`

## Goal

After `B3`, `B4`, `N118`, `N163`, and `N234`, the strongest honest reopening
route is no longer:
- observer-side global promotion,
- or silent promotion of downstream asymmetry into selector closure.

The strongest honest reopening route is narrower:
- future derivation of a non-trivial source-topology invariant of the kernel or
  source limit,
- followed by a strict-core, basis-independent promotion of that invariant into
  a selector datum upstream of the observer.

`T14` does not claim that such a theorem is already proved.

`T14` does something narrower:
- writes a packet-ready future theorem spec for a `Source Topology Selector`
  route,
- isolates the minimal assumptions under which this route would make sense,
- makes explicit that the observer may appear only as a downstream algebraic
  pushforward witness,
- preserves the no-false-pass boundary around `QW-2191`.

## Target theorem

### Informal statement

If the strict core exports a non-trivial source-topology invariant at the
kernel or source limit, and if that invariant admits a quotient-safe,
basis-independent promotion into a selector datum before observer stages, then
that invariant may serve as a `Source Topology Selector` for the declared
scope.

In that case the downstream observer asymmetry may be read only as an algebraic
pushforward witness of the upstream selector source, not as the cause of the
selector.

### Formal target

```text
T14_SourceTopologySelector_Theorem

Assume:
  A1. a strict-core exported non-trivial source-topology invariant
      tau_src exists at the declared kernel/source limit;
  A2. tau_src is derived from declared strict-core source objects only and does
      not import an external selector anchor;
  A3. a basis-independent strict-core promotion map
      Pi_sel : tau_src -> selector_datum
      is exported;
  A4. Pi_sel is quotient-safe with respect to the present O(2) uniqueness
      frontier from QW-2191 and does not rely on observer-side promotion;
  A5. the promoted selector datum is exported upstream of the observer and is
      compatible with the declared order
      nadsoliton -> light -> matter -> emergent observer;
  A6. downstream observer asymmetry is used only as algebraic pushforward
      witness and not as primary selector source;
  A7. all claims remain inside the anti-overclaim boundary.

Then:
  tau_src may serve as a strict-core Source Topology Selector in the declared
  scope, and the promoted selector datum may be used as the upstream witness
  resolving the current selector-source gap for that scope.
```

## Strict-admissible support

1. `QW-2191`
   - strict-core uniqueness obstruction, O(2) family, no kernel-alone
     uniqueness
2. `QW-2206`
   - local topological protection layer
3. `B3`
   - topological selector bridge packet
4. `B4`
   - minimal `sigma_int_candidate`
5. `N118`
   - selector or symmetry-breaking requirement remains active on current repo
     state
6. `N163`
   - observer information deficit is downstream symptom, not primary source
7. `N234`
   - exported emergent-observer chain does not justify global selector closure
     or global `QW-2191` discharge
8. `A10`
   - anti-overclaim boundary

## Minimal lemma DAG

### L1. A non-trivial source-topology invariant exists

```text
L1:
The strict core exports a non-trivial topological invariant tau_src at the
declared kernel/source limit.
```

Support now:
- not discharged

Motivation:
- `QW-2206`
- `B3`
- `B4`

### L2. The invariant is upstream of the observer

```text
L2:
tau_src is exported on the source side and remains upstream of all
observer/coarse-graining/readout constructions.
```

Support now:
- only route discipline

Motivation:
- `N163`

### L3. The invariant admits a basis-independent selector promotion

```text
L3:
There exists a strict-core map Pi_sel taking tau_src to a selector datum
without dependence on an arbitrary basis choice.
```

Support now:
- not discharged

Motivation:
- `B3`

### L4. The promotion is quotient-safe at the QW-2191 frontier

```text
L4:
The promoted selector datum resolves the declared O(2) ambiguity in a
quotient-safe way and therefore may count as a selector witness for the
declared scope.
```

Support now:
- not discharged

Motivation:
- `QW-2191`
- `N118`

### L5. Observer asymmetry is only a downstream pushforward witness

```text
L5:
Any downstream observer asymmetry may be used only as an algebraic pushforward
of the upstream selector datum and cannot serve as the primary selector source.
```

Support now:
- supported as a boundary statement

Motivation:
- `N163`
- `N234`

### L6. L1-L5 yield a Source Topology Selector theorem

```text
L6:
If L1-L5 hold, then tau_src can serve as a Source Topology Selector in the
declared strict-core scope.
```

Support now:
- theorem-spec only

## Minimal assumption map

### Source-side assumptions

1. a non-trivial source-topology invariant is actually exported;
2. the invariant is strict-core only and not an external selector import;
3. the invariant is not merely a control-lane fit or observer-side summary.

### Promotion assumptions

1. a selector promotion map exists;
2. the promotion is basis-independent;
3. the promotion is quotient-safe across the `QW-2191` ambiguity family;
4. the promoted datum is exported before observer stages.

### Observer-side assumptions

1. downstream asymmetry remains a witness only;
2. observer-side stability cannot be used as the missing upstream source by
   itself;
3. no global discharge is claimed from observer-side closure candidates alone.

## Acceptance skeleton

The theorem spec is acceptable only if all of the following stay explicit:

1. the theorem is future-route only and not a current discharge;
2. a non-trivial topological invariant must be explicitly exported, not only
   suggested heuristically;
3. basis-independence must be proved, not inferred from one successful chain;
4. observer-side asymmetry remains downstream witness only;
5. no step is described as current global selector closure;
6. no step is described as current global `QW-2191` discharge;
7. no step silently treats `sigma_int_candidate` as already sufficient.

## What this theorem would establish if discharged

If discharged, `T14` would establish:
- a strict-core upstream selector source based on source topology;
- a basis-independent promotion route from source topology to selector datum;
- a route to resolving the `QW-2191` uniqueness frontier in the declared
  scope;
- a principled reason why downstream observer asymmetry is evidential rather
  than causal.

It would not by itself establish:
- legacy-to-strict kernel bridge;
- full ToE closure;
- every future global extension of the theory;
- unrestricted global discharge beyond the declared scope.

## Residual blockers after the spec

Even after `T14` is written, the following remain open:

- actual export of a non-trivial `tau_src`;
- actual basis-independent promotion map `Pi_sel`;
- actual quotient-safe proof at the `QW-2191` frontier;
- actual theorem-level integration with the existing positive preobserver chain;
- any current global discharge statement.

## Anti-overclaim

`T14` does not claim:
- theorem-level PASS,
- current selector closure,
- current global selector closure,
- current `QW-2191` discharge,
- current proof that the exported observer lane already globalizes,
- current proof that `sigma_int_candidate` already solves the selector problem.

## Product of the step

- packet-ready future theorem spec for a `Source Topology Selector` route;
- explicit separation of upstream source-topology promotion from downstream
  observer witness;
- explicit no-false-pass barrier against premature global promotion.

## Next step

Natural next move:
- attempt a future route packet that exports one candidate
  `tau_src` at the source/kernel limit,
- but keep it explicitly below theorem-level discharge until
  basis-independence and `QW-2191` quotient-safety are actually shown.
