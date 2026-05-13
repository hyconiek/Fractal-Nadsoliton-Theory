# P1417 PC8 Noncyclic Selector-Margin Amplifier Design Packet (EN/PL)

Status: `P1417_EXECUTED_PC8_DESIGN_FREEZE_NO_FALSE_PASS`
As of: `2026-05-13`

## Professorial decision

After `P1416 FAIL_PC7_MARGIN`, move to a noncyclic `PC8` provider with explicit margin-amplifier design.

## Execution

- Script: `p1417_pc8_noncyclic_selector_margin_amplifier_design_checkpoint.py`
- Artifact: `generated/p1417_pc8_noncyclic_selector_margin_amplifier_design_summary.json`

## Result

Design frozen for `PC8` with stricter selector-margin target (`>= 0.0024`) while keeping
transport and replay thresholds unchanged.

No bridge claim, no promotion claim.

## Lay explanation (PL)

Po ludzku: poprzednia wersja była za „mało wyraźna”, więc projektujemy nową,
która ma mocniej odróżniać najlepszy wybór od reszty, ale nadal musi być stabilna.

## Recommendation

Execute `P1418` (first PC8 transport + replay run) and accept only threshold-based verdict.
