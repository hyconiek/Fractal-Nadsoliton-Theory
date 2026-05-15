# P1806 S756 Strict Unified Nonproxy EA-EH-ELg Residual Runpack Contract Packet (PL)

Status: `P1806_EXECUTED_STRICT_UNIFIED_NONPROXY_EA_EH_ELG_RESIDUAL_RUNPACK_CONTRACT_PACKET_NO_FALSE_PASS`

## Cel

Wykonać następny uczciwy krok po `P1805`:

zdefiniować jeden wspólny run-pack kontrakt dla
`E_A^mu + E_H + EL_g` na tym samym freeze/manifeście,
aby binarnie oceniać wejście do `TG1_BW`.

To NIE jest claim PASS. To jest warstwa wykonawcza + walidacyjna.

## Kontrakt run-pack (minimum)

Run-pack jest `ADMISSIBLE` tylko jeśli zawiera jednocześnie:

1. `shared_background_family_id` (identyczne dla EA/EH/ELg),
2. `boundary_control_clause_id` (powiązane z `P1761/P1762`),
3. `weak_form_h1_4d_ready = true`,
4. trzy sekcje residual:
   - `ea_mu_residual`,
   - `eh_residual`,
   - `elg_residual`,
5. jawny divergence trace dla każdego residualu, gdy nie jest zero,
6. jawny werdykt lokalny per residual: `PASS_ZERO` albo `OPEN_OBSTRUCTION_WITH_TRACE`.

## Reguła binarna TG1

`TG1_BW = PASS_ZERO` tylko jeśli:

- `ea_mu_residual.verdict = PASS_ZERO`,
- `eh_residual.verdict = PASS_ZERO`,
- `elg_residual.verdict = PASS_ZERO`,
- wszystkie trzy na tym samym `shared_background_family_id` i `manifest_id`.

W przeciwnym razie:

- `TG1_BW = OPEN_LOCKED_BY_UNIFIED_NONPROXY_RESIDUAL`.

## Explicit anti-false-pass guard

Niedozwolone:

1. PASS dla dwóch residuali + brak trzeciego.
2. PASS sektorowy (`EA/EH`) podniesiony do BW global.
3. Mieszanie różnych freeze/background family między residualami.

## Co zostało dowiedzione

1. Mamy teraz formalny, maszynowo walidowalny kontrakt jednego run-packu.
2. Kryterium odblokowania BW jest jednoznaczne i binarne.
3. Lock BRST/CUT pozostaje poprawnie wymuszany przy `TG1 != PASS_ZERO`.

## Co pozostaje OPEN

1. Rzeczywisty kompletny run-pack z trzema `PASS_ZERO`.
2. Theorem-level witness dla BW.
3. Dalsze `TG2_BRST` i `TG3_CUT` po BW.

## Produkt

- Packet + validator + template JSON + checkpoint JSON.
