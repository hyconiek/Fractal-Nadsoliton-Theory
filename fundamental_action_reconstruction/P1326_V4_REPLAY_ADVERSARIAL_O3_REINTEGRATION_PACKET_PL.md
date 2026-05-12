# P1326 — V4 replay/adversarial + O3 reintegration packet (PL)

Status: `V4_REINTEGRATION_READY_THEOREM_PENDING`
As of: `2026-05-12`
Depends on: `P1325`

## Cel
Wykonać następny uczciwy krok: niezależny replay i adversarial dla `v4`,
a następnie ocenić gotowość do reintegracji z formalnym O3-EXCLUSION.

## Artefakt wykonawczy
- skrypt: `p1326_v4_replay_adversarial_reintegration.py`
- raport: `generated/p1326_v4_replay_adversarial_reintegration_report_v1.json`

## Wynik
- replay (3 seedy): `PASS` (we wszystkich runach `sign_diversity=2`, `flip_rate<=0.08`),
- adversarial edge attack: `PASS` (`unstable_rate=0.0128`),
- `o3_reintegration_ready = true`.

## Decyzja profesorska
To jest istotny postęp: kandydat `v4` przeszedł niezależną walidację
numeryczno-adversarialną i może zostać formalnie wpięty do O3.

Jednocześnie: nadal brakuje theorem-level discharge; dlatego:
- `strict_core_selector_source_exported = false`,
- `QW-2191` strict nadal `NOT_CLOSED`.

## Konsekwencja
Program przeszedł z etapu "kandydat testowy" do etapu
"kandydat gotowy do formalizacji twierdzenia".

## Czego dokument nie twierdzi
- Nie twierdzi strict closure.
- Nie twierdzi ToE closure.
- Nie zastępuje dowodu formalnego samą walidacją numeryczną.

## Rekomendowany następny uczciwy krok
Uruchomić **P1327 theorem-level O3-EXCLUSION update with v4**:
1. sformalizować lemma `V4_selector_consistency`,
2. wykazać eliminację residual loopholes w logice O3,
3. dopiero wtedy rozstrzygać `QW-2191` strict status.

## Dla laika
Nowy kompas `v4` przeszedł powtórne i "złośliwe" testy bardzo dobrze.
To znaczy, że jest mocnym kandydatem. Ale w nauce to jeszcze pół kroku:
trzeba go teraz zamknąć formalnym dowodem, zanim ogłosimy pełne domknięcie.
