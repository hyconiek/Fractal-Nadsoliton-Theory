# P1625 / S575 — Previous candidate failure audit + next strict move

## Cel
Sprawdzić poprzednie badania i potwierdzić, które kandydaty strict-core nie domknęły się,
aby nie powtarzać tych samych prób i przejść do kolejnego uczciwego ruchu.

## Wejścia
- `generated/p1583_s533_formal_proof_object_and_global_stability_composition_summary.json`
- `generated/p1585_s535_formal_l1_l3_theorem_objects_and_toe_composition_summary.json`
- `generated/p1596_s546_selector_uniqueness_bridge_theorem_object_summary.json`
- `generated/p1605_s555_np1_provider_instantiation_and_replay_summary.json`
- `generated/p1623_s573_strict_selector_uniqueness_theorem_object_and_variational_log_summary.json`
- `generated/p1624_s574_noncyclic_selector_witness_from_strict_kernel_summary.json`

## Wyjście
- `generated/p1625_s575_previous_candidate_failure_audit_and_next_strict_move_summary.json`

## Rygor
- strict-only, bez legacy bridge.
- bez walidacji zewnętrznej.
- full chain: `K_strict -> coeff -> L_total -> EOM` jako jedyna baza dalszych decyzji.
