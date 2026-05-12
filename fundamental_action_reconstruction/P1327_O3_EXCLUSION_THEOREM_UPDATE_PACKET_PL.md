# P1327 — O3-EXCLUSION theorem update packet with v4 (PL)

Status: `THEOREM_OBLIGATIONS_PARTIAL_NUMERICS_PASS_FORMAL_EXPORT_PENDING`
As of: `2026-05-12`
Depends on: `P1326`

## Cel
Wykonać krok formalizacji po `P1326`: przejść od walidacji numerycznej `v4` do
jawnej tabeli obligations theorem-level dla O3-EXCLUSION.

## Artefakt wykonawczy
- skrypt: `p1327_o3_exclusion_theorem_obligation_checker.py`
- raport: `generated/p1327_o3_exclusion_theorem_obligation_checker_report_v1.json`

## Wynik
- obligations numeryczne: `PASS` (`V4_replay_pass`, `V4_adversarial_pass`,
  `V4_reintegration_ready`),
- obligations formalne: `PENDING`
  (`theorem_V4_selector_consistency_formally_exported = false`,
  `residual_loophole_elimination_formally_exported = false`),
- `theorem_ready = false`,
- `qw2191_strict_status = NOT_CLOSED`.

## Decyzja profesorska
W tej chwili mamy mocny materiał empiryczny, ale brakuje dwóch kluczowych
obowiązków dowodowych theorem-level. To oznacza, że strict closure jest nadal
niedopuszczalne.

## Konsekwencja
- Bramka formalna O3 nie jest jeszcze domknięta.
- Następny ruch musi być czysto dowodowy, nie kolejny tuning numeryczny.

## Czego dokument nie twierdzi
- Nie twierdzi strict closure.
- Nie twierdzi ToE closure.
- Nie twierdzi, że numerics wystarczą bez formal export.

## Rekomendowany następny uczciwy krok
Uruchomić **P1328 formal export packet** obejmujący:
1. lemma `V4_selector_consistency` (pełny statement + proof sketch),
2. lemma `Residual_loophole_elimination` (jawna eliminacja `open(Z2/eps)`),
3. aktualizację decyzji O3 na bazie tych dwóch formalnych eksportów.

## Dla laika
Przetestowaliśmy kompas bardzo dobrze i działa stabilnie. Ale to jeszcze nie
koniec, bo w nauce potrzebny jest też formalny dowód "dlaczego". Mamy więc
solidne dane, ale brakuje dwóch ostatnich cegiełek dowodowych.
