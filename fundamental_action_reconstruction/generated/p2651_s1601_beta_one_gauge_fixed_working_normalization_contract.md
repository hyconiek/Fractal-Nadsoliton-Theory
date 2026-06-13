# P2651/S1601 beta=1 gauge-fixed working-normalization contract

Status: `P2651_BETA_ONE_GAUGE_FIXED_WORKING_NORMALIZATION_CONTRACT_NO_SOURCE_NO_FALSE_PASS`

## Content-first anti-duplication audit

This packet greps gauge-fixed beta, rescaling invariants, role boundaries, and source nonclosure before adding the contract.

- `gauge_fixed_beta_content`: 54 hits
- `rescaling_invariant_content`: 3652 hits
- `role_boundary_content`: 1604 hits
- `source_nonclosure_content`: 13900 hits

## Orbit invariance

d' = a*d, beta' = beta/a^eta preserves x=beta*d^eta and the envelope 1/(1+x).
For beta>0 choose a=beta^(1/eta), giving beta'=1 and d'=beta^(1/eta)*d.
Max envelope error on audited grid: `1.388e-17`.

## Tail-ratio gauge warning

Tail ratios are gauge-respecting only when the distance coordinates are scaled with beta; setting beta=1 while leaving raw distances unchanged is a new gauge choice, not an invariant operation.
Max gauge-respecting ratio error: `2.776e-17`.
Max raw-distance substitution error: `0.699767037772`.

## Claim boundary matrix

| claim | allowed? | status | required label |
| --- | ---: | --- | --- |
| `strict_beta_equals_one_working_gauge` | `True` | `ALLOWED_AS_GAUGE_FIXED_WORKING_NORMALIZATION` | must say beta=1 is a chosen representative of the damping orbit |
| `target_independent_beta_source` | `False` | `BLOCKED_BY_P2649_P2650` | needs canonical length/UV unit plus source identity |
| `modified_compressed_successor_semantics` | `True` | `ALLOWED_AS_GAUGE_DECLARED_DESCRIPTION` | must declare the beta=1 gauge and raw distance convention |
| `unchanged_legacy_inverse_hierarchy_transfer` | `False` | `REJECTED` | P2643-P2645 rejection remains active |
| `blind_holdout_empirical_support` | `True` | `ALLOWED_ONLY_IF_DATA_GAUGE_AND_UNIT_MAP_ARE_DECLARED` | cannot replace beta source theorem |
| `role_bearing_ltotal_damping_term` | `False` | `BLOCKED` | needs beta source, typed unit, role-transfer, and selector/source discharge |

## Verdict

P2651 is the honest fallback after P2649-P2650: beta=1 may be used as a declared gauge-fixed working normalization, because every positive beta lies on a rescaling orbit with a beta=1 representative.  The contract prevents that representative from being promoted to a beta source, legacy role transfer, or role-bearing dynamics.

Decision: `BETA_ONE_GAUGE_FIXED_WORKING_NORMALIZATION_CONTRACT_ONLY__NO_SOURCE_NO_LTOTAL_NO_TOE`.
Beta=1 working gauge allowed? `True`.
Beta source exported now? `False`.
Role-bearing L_total now? `False`.
ToE closure now? `False`.

## Professorial closure path

- Use beta=1 only with an explicit gauge declaration and transformed distance convention.
- For empirical P2647/P2648 runs, require a data-unit map into the declared beta=1 gauge before comparing tail ratios.
- For theory closure, the missing next atom is still a typed nadsoliton metric/UV source, not a larger numerical scan.
- Rerun role-transfer only after beta source, bridge completion, and selector/source blockers are independently discharged.

## Next honest step

Adopt the P2651 gauge-fixed contract in downstream packets, then either build a typed nadsoliton metric/UV source theorem or execute P2647/P2648 with a declared data-unit map as empirical support only.

## Negative exports

- `canonical_length_uv_unit_exported`: `False`
- `target_independent_beta_identity_exported`: `False`
- `beta_source_exported`: `False`
- `legacy_role_transfer_revalidated`: `False`
- `empirical_holdout_promoted_to_source`: `False`
- `bridge_completion_exported`: `False`
- `q_w_2191_discharged`: `False`
- `role_bearing_ltotal_reenabled`: `False`
- `toe_closure_claimed`: `False`
