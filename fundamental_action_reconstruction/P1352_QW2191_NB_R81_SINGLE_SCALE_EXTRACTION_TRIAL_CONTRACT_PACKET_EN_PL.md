# P1352 QW-2191 NB R8.1 Single-Scale Extraction Trial Contract Packet (EN+PL)

Status: `P1352_CONTRACT_FROZEN_NO_FALSE_PASS`
As of: `2026-05-12`
Depends on: `P1347`, `P1348`, `P1349`, `P1351`
Artifact: `generated/p1352_qw2191_nb_r81_single_scale_extraction_trial_contract_summary.json`

## Goal

Move from documentation-only intent to an executable strict contract object for:

- NB-lane QW-2191 closure discipline,
- single-scale extraction bridge to SM/GR observables,
- mandatory external blind-audit governance.

## Executed step

`p1352_qw2191_nb_r81_single_scale_extraction_trial_contract_checkpoint.py`
was executed and exported one machine-readable trial contract summary.

## Frozen contract (strict scope)

1. One declared scale (`mu_*`) and one declared renormalization scheme.
2. Target observable set: `(g1, g2, g3)` + one gravity effective observable.
3. Mandatory outputs: public residual table, PASS/FAIL status, incident log.
4. Governance hard flags:
   - no new axioms,
   - no silent legacy-role transfer,
   - blind audit required,
   - rollback on fail required.

## Professor decision

The next honest move is not another narrative packet, but direct execution of
`R8_1_BLIND_AUDIT_EXECUTION_AND_RESIDUAL_PUBLICATION` on the frozen contract.

## PL — dla laika

Zrobiliśmy krok „z planu do działania”:

- nie tylko mówimy co policzyć,
- ale zapisaliśmy to jako konkretny kontrakt techniczny,
- który ma dać jasny wynik: zgodność albo niezgodność z danymi.

To jest uczciwy naukowo test: jeśli nie przejdzie, teoria wraca do poprawy.
