# P1802 S752 Strict Cutkosky Unitarity Witness Intake And TG3 Gate Link Packet (PL)

Status: `P1802_EXECUTED_STRICT_CUTKOSKY_UNITARITY_WITNESS_INTAKE_AND_TG3_GATE_LINK_PACKET_NO_FALSE_PASS`
As of: `2026-05-15`

## Cel

Domknąć ostatni formalny intake-link na ścieżce `BW -> BRST -> Cutkosky`:

1. zdefiniować minimalny witness intake dla unitarności cięć,
2. związać intake z `TG3_CUTKOSKY_GLOBAL_UNITARITY`,
3. zablokować claim TG3 PASS bez jawnych discontinuity/ghost-pole witnessów.

## Zakres strict-only

Brak legacy bridge i brak shortcutów. To wyłącznie intake+gate-link dla TG3.

## Minimalny witness pack TG3

Wymagane artefakty:

1. `cut_discontinuity_witness_id` (jawne cięcia amplitud),
2. `optical_theorem_match_id` (zgodność części urojonej i sumy cięć),
3. `ghost_pole_exclusion_id` (brak niedozwolonych biegunów ghost),
4. `background_family_consistency_id`,
5. `shared_freeze_id` zgodny z BW+BRST,
6. `check_digests` dla wszystkich punktów.

## TG3 gate-link policy

`TG3_CUTKOSKY_GLOBAL_UNITARITY = PASS` tylko jeśli:

- `G_BW = PASS_ZERO`,
- `TG2_BRST_GLOBAL_NILPOTENCY = PASS`,
- witness pack TG3 kompletny,
- optical-theorem match jest jawnie spełniony,
- ghost-pole exclusion jest jawnie spełnione,
- brak freeze mismatch.

W przeciwnym razie:

- `TG3_CUTKOSKY_GLOBAL_UNITARITY = OPEN_LOCKED_BY_UNITARITY_WITNESS`.

## Co jest dowiedzione

1. TG3 ma pełny format intake dowodowego zamiast statusu deklaratywnego.
2. Kolejność `BW -> TG2 -> TG3` jest formalnie egzekwowalna.
3. Nie można ogłosić unitarności bez jawnych witnessów discontinuity/optical/ghost.

## Co pozostaje OPEN

1. Realny TG3 witness pack na wspólnym freeze.
2. Faktyczna theorem-level unitarity closure dla miksu `spin-2 + SM`.

## Produkt

- Packet PL.
- Checkpoint JSON policy + TG3 witness intake template.
