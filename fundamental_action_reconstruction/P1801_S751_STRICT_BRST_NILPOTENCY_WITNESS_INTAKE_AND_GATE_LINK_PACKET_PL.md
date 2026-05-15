# P1801 S751 Strict BRST Nilpotency Witness Intake And Gate Link Packet (PL)

Status: `P1801_EXECUTED_STRICT_BRST_NILPOTENCY_WITNESS_INTAKE_AND_GATE_LINK_PACKET_NO_FALSE_PASS`
As of: `2026-05-15`

## Cel

Zamknąć lukę między `BW PASS` a `TG2_BRST_GLOBAL_NILPOTENCY`:

1. zdefiniować minimalny witness intake dla BRST nilpotency,
2. jawnie powiązać ten intake z gate lock matrix (`TG2`, `SV7`),
3. zablokować claim BRST PASS bez jawnego `Q^2=0` witness pack.

## Zakres strict-only

Brak legacy bridge i brak shortcutów. To tylko warstwa intake + gate-link.

## Minimalny witness pack BRST

Wymagane artefakty:

1. `brst_charge_definition_id` (jawny operator `Q`),
2. `nilpotency_check_id` (jawny check `Q^2`),
3. `cohomology_subspace_check_id` (fizyczna podprzestrzeń),
4. `ghost_sector_consistency_id`,
5. `shared_freeze_id` zgodny z BW run,
6. `check_digests` dla wszystkich powyższych.

## TG2 gate-link policy

`TG2_BRST_GLOBAL_NILPOTENCY = PASS` tylko jeśli:

- `G_BW = PASS_ZERO`,
- witness pack BRST kompletny,
- `Q^2` jawnie uproszczone do zera,
- brak freeze mismatch.

W przeciwnym razie:

- `TG2_BRST_GLOBAL_NILPOTENCY = OPEN_LOCKED_BY_NILPOTENCY_WITNESS`.

## Co jest dowiedzione

1. BRST gate ma jawny format dowodu zamiast deklaratywnego statusu.
2. Zależność `BW -> TG2` jest technicznie egzekwowalna.
3. Brak możliwości "paper PASS" dla BRST bez `Q^2=0` evidence.

## Co pozostaje OPEN

1. Realny nonproxy BRST witness pack na wspólnym freeze z BW.
2. Pozytywny werdykt TG2 i dopiero potem wejście w Cutkosky TG3.

## Produkt

- Packet PL.
- Checkpoint JSON policy + BRST witness intake template.
