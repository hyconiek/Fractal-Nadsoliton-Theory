# P1795 S745 Strict Nonproxy Covariant Export + BW/BRST/CUT Gate Readiness Packet (PL)

Status: `P1795_EXECUTED_STRICT_NONPROXY_COVARIANT_EXPORT_BW_BRST_CUT_GATE_READINESS_PACKET_NO_FALSE_PASS`
As of: `2026-05-15`

## Cel

Zamknąć "następny uczciwy krok" po `P1790..P1794` w trybie strict-only:

1. zebrać obecne świadectwa dla `E_Aμ`, `E_H`, `EL_g`,
2. odseparować `LOCAL` od `GLOBAL` oraz `SCAFFOLD` od `FULL EXPORT`,
3. wyprowadzić minimalny, niesprzeczny kontrakt wejściowy dla bramki `G_BW`,
4. utrzymać blokadę `G_BRST` i `G_CUT`, dopóki `G_BW != PASS_ZERO`.

Pakiet nie claimuje theorem-level closure.

## Zakres strict-only

Ten pakiet używa wyłącznie strict-side artefaktów i nie wprowadza żadnego
legacy->strict bridge. Nie przenosi także legacy physical-role claims.

## Aktualna mapa dowodowa (stan uczciwy)

### 1) Covariant nonproxy `E_Aμ`

- Istnieje jawny eksport kandydacki (`P1764`) oraz ścieżka dostaw `W1`.
- Status: `LOCAL_EXPORT_PRESENT`, `GLOBAL_GATE_NOT_CLOSED`.

### 2) Covariant nonproxy `E_H`

- Istnieje jawny eksport kandydacki (`P1764`) + boundary-control clauses (`P1761/P1762`).
- Status: `LOCAL_EXPORT_PRESENT`, `BOUNDARY_SAFETY_PARTIAL`,
  `GLOBAL_GATE_NOT_CLOSED`.

### 3) Metric `EL_g`

- Istnieje jawny eksport kandydacki (`P1765`) oraz komponentowy scaffold (`P1718..P1727`).
- Status: `LOCAL_EXPORT_PRESENT`, `FULL_TENSOR_CLOSURE_OPEN`.

### 4) Bianchi/Ward

- Sekwencja i lock bramkowy są jawnie zdefiniowane (`P1767`, `P1788`).
- Status: `INPUT_CONTRACT_READY`, `PASS_WITNESS_OPEN`.

### 5) BRST/Cutkosky

- Obie bramki pozostają downstream i zablokowane bez `G_BW: PASS_ZERO`.
- Status: `LOCKED_BY_POLICY`.

## Minimalny kontrakt wejściowy do uczciwego uruchomienia `G_BW`

Wymagane jednocześnie (ten sam freeze i ta sama rodzina tła):

1. `EA_covariant_nonproxy` dostarczone jako expression pack + metadata,
2. `EH_covariant_nonproxy` dostarczone jako expression pack + boundary clause id,
3. `ELg_nonproxy` dostarczone jako expression pack + konwencje indeksowe,
4. boundary-control witness (`BC_family_id`, `falloff_clause`, `boundary_mode`),
5. intake walidowany przez `P1793` bez freeze mismatch.

Jeśli którykolwiek punkt niespełniony -> verdict: `OPEN_OBSTRUCTION_WITH_TRACE`.

## No-false-pass werdykt

Dopuszczalne tylko dwa wyniki runtime:

1. `PASS_ZERO` (jawny residual/witness),
2. `OPEN_OBSTRUCTION_WITH_TRACE` (jawny blocker + trace).

Niedopuszczalne:

- PASS bez residualu,
- PASS lokalny podniesiony do globalnego,
- BRST/Cutkosky start przed `G_BW: PASS_ZERO`.

## Update reguły stanu

Po każdym runie:

1. walidacja intake (`P1793`),
2. update wektora stanu (`P1794`),
3. dopiero wtedy ocena, czy `G_BW` zmienia się na `PASS_ZERO`.

`G_BRST` i `G_CUT` pozostają `LOCKED`, jeżeli `G_BW != PASS_ZERO`.

## Co jest dowiedzione w tym pakiecie

1. Kontrakt wykonawczy dla uczciwego `G_BW` jest kompletny proceduralnie,
2. rozróżnienie `LOCAL/GLOBAL` i `SCAFFOLD/FULL` jest jawnie utrzymane,
3. nie ma logicznej luki pozwalającej na "silent gate promotion".

## Co pozostaje open

1. rzeczywisty witness `PASS_ZERO` dla `G_BW`,
2. pełna tensorowa domkniętość metryczna,
3. theorem-level BRST nilpotency,
4. theorem-level Cutkosky unitarity,
5. renormalization/counterterm closure jako globalny witness.

## Produkt

- Ten packet PL.
- Checkpoint JSON do automatycznej konsumpcji i dalszego update stanu.
