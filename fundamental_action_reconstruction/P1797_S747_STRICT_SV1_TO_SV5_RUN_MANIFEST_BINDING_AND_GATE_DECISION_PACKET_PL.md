# P1797 S747 Strict SV1..SV5 Run Manifest Binding And Gate Decision Packet (PL)

Status: `P1797_EXECUTED_STRICT_SV1_TO_SV5_RUN_MANIFEST_BINDING_AND_GATE_DECISION_PACKET_NO_FALSE_PASS`
As of: `2026-05-15`

## Cel

Dostarczyć brakujący element wykonawczy między template ledgera (`P1796`) a werdyktem bramkowym:

1. jawne powiązanie `SV1..SV5` z artefaktami `EA/EH/ELg/boundary/H1`,
2. deterministyczna decyzja `G_BW` na podstawie kompletności i witnessów,
3. twarde utrzymanie locków `G_BRST/G_CUT` gdy `G_BW != PASS_ZERO`.

## Zakres strict-only

Pakiet używa wyłącznie strict-side artifact binding. Brak legacy bridge, brak proxy,
brak global claim bez theorem witness.

## Manifest binding (wymagane pola)

Manifest runu musi zawierać:

- `freeze_id_common`,
- `sv_bindings`:
  - `SV1 -> EA_covariant_nonproxy_artifact_id`,
  - `SV2 -> EH_covariant_nonproxy_artifact_id`,
  - `SV3 -> ELg_nonproxy_artifact_id`,
  - `SV4 -> boundary_control_witness_id`,
  - `SV5 -> H1_4D_weak_form_witness_id`,
- `validation`:
  - `intake_validation_result_id` (z `P1793`),
  - `state_update_protocol_id` (z `P1794`),
- `gate_decision_inputs`:
  - `bw_residual_report_id`,
  - `bw_check_digest`.

## Gate decision policy

1. `G_BW = PASS_ZERO` tylko jeśli:
   - wszystkie bindingi obecne,
   - ten sam `freeze_id_common`,
   - `P1793` = PASS,
   - `bw_residual_report` jawnie zero + digest.
2. W każdym innym przypadku:
   - `G_BW = OPEN_OBSTRUCTION_WITH_TRACE`.
3. Gdy `G_BW != PASS_ZERO`:
   - `G_BRST = LOCKED`,
   - `G_CUT = LOCKED`.

## Co jest dowiedzione

1. Decyzja `G_BW` ma teraz explicit machine policy bez szarej strefy.
2. Brak możliwości "przeskoku" do BRST/CUT przez niekompletny manifest.
3. Łańcuch `P1793 -> P1794 -> G_BW decision` jest egzekwowalny.

## Co pozostaje OPEN

1. Realny manifest z pełnymi artifact IDs i zerowym raportem BW.
2. Theorem-level BRST nilpotency witness.
3. Theorem-level Cutkosky unitarity witness.
4. Global Helmholtz integrability closure.

## Produkt

- Packet PL.
- Checkpoint JSON policy + runnable manifest template JSON.
