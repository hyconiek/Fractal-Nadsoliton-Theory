# Strict kernel legacy role-transfer audit draft

Status: `DRAFT_V0_PRE_AUDIT__ALL_LEGACY_ROLES_BLOCKED_PENDING_STRICT_THEOREMS`

This file starts the release scaffold for updating legacy physical-role claims after the legacy -> strict bridge. The purpose is not to preserve legacy roles by inertia. The strict kernel is richer than the legacy kernel, so each role must be re-audited against strict nonlinear compression, phase/topological data, and selector/source constraints.

## 1. Audit principle

Legacy roles may be transferred to `K_strict_gate` only after both conditions hold:

1. the bridge `K_legacy_ont -> K_strict_gate` is theorem-grade rather than finite-comparison scope;
2. a separate role-transfer theorem proves which physical meaning survives unchanged, which is modified, and which is rejected.

Until then, the default verdict is blocked.

## 2. Role-transfer table

| Legacy role | Legacy formula / candidate | Strict-side update requirement | Current verdict |
|---|---|---|---|
| Electroweak angle role | `sin^2(theta_W)=alpha_geo/12` | Needs an `alpha_geo` electroweak-role theorem after strict amplitude normalization; scalar normalization alone is not a physical-role proof. | blocked |
| Electromagnetic inverse role | `alpha_EM^-1 = alpha_geo/(2*beta_tors)*(1-beta_tors)` | Needs both an `alpha_geo` role theorem and a `beta_tors` strict-role theorem; strict damping has `beta=1`, `eta=9/5`, so the legacy linear torsion denominator cannot be copied. | blocked / likely modified if it survives |
| Gravity hierarchy role | `beta^N` hierarchy | Needs a strict beta-power hierarchy successor theorem; nonlinear compression `d^eta` changes the semantics of a pure beta-power legacy hierarchy. | blocked / likely modified if it survives |
| Selector candidate | `beta_tors -> chi_11` | Needs a theorem reducing continuous legacy torsion to the GF(2) phase/topology selector source; current audits identify `chi_11` as a bottleneck but not as derived. | blocked |

## 3. Strict characteristics that force role updates

The following strict characteristics are not optional decorations. They can change physical-role semantics:

- nonlinear damping/compression `1 + beta*d^eta` with audited `beta=1`, `eta=9/5`;
- GF(2) phase-sign solution with four forced flip edges;
- path/cycle cohomology distinction between `C0` anchor and closed-cycle `C1/im(delta)` generator;
- finite carrier/node/edge bookkeeping;
- open selector/source atom for `chi_11`;
- open strict source atoms for APD, phase/frequency, and damping parameters.

## 4. Minimal obligation atoms

The current minimal obligation lattice records four missing theorem atoms before legacy role transfer can be claimed:

```text
alpha_geo electroweak role
beta_tors strict role
beta-power hierarchy successor
chi11 selector/source
```

The pre-audit matrix blocks all four audited roles under the current all-open frontier. The theorem-frontier reports further show that low-weight additions can unlock selector-only signatures, but no singleton or pair closes bridge, role transfer, or ToE.

## 5. Updated strict-role target semantics

The next release should not ask only, "does the legacy formula survive?" It should ask:

1. **unchanged transfer:** does the exact legacy expression remain valid after strict completion?
2. **modified transfer:** does a strict successor expression exist after adding nonlinear compression/topology?
3. **rejection:** was the legacy role merely an artifact of the incomplete kernel?

For example, if the electromagnetic inverse role survives, it likely cannot remain a raw function of `beta_tors` alone unless a theorem explains how `beta_tors` maps into the strict `beta/eta` compression and the `chi_11` phase selector. If the gravity hierarchy survives, it likely needs a strict hierarchy using compression data rather than silent `beta^N` inheritance.

## 6. Hard limits

- No legacy role is transferred in this draft.
- No `alpha_geo` electroweak role theorem is claimed.
- No strict `beta_tors` role theorem is claimed.
- No `beta_tors -> chi_11` theorem is claimed.
- No strict beta-power hierarchy successor is claimed.
- No QW-2191 selector discharge is claimed.
- No ToE closure is claimed.

## 7. Next release-build tasks

1. Turn each role row into a theorem target with assumptions and strict outputs.
2. Distinguish unchanged survival, strict modification, and rejection for each role.
3. Add explicit dependency links from role-transfer targets to bridge-source atoms and `chi_11` selector-source atoms.
4. Only after those theorem targets pass, update the release tex/md physical-role narrative.
