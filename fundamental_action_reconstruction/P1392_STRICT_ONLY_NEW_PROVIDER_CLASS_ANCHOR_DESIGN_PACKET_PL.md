# P1392 Strict-only New Provider-Class / Anchor Design Packet (No Legacy Bridge) — PL

Status: `P1392_EXECUTED_NEW_PROVIDER_CLASS_DESIGN_BASELINE_FROZEN`
As of: `2026-05-13`

## Cel

Po formalnym obstruction export (`P1391`) zaprojektować nową klasę provider/anchor,
aby otworzyć niecykliczną gałąź dalszych prób dla
`F_Nadsoliton => L_SM + L_GR` bez mostu do legacy.

## Decyzja profesorska (noncyclic)

Zamykamy starą klasę lokalnych patch-weight loopów (`ce6` obstruction v1)
i startujemy nową klasę providerów `PC2` z odseparowanym baseline.

## Definicja `PC2` (v1)

- `provider_class_id = PC2_strict_boundary_anchor_v1`
- `anchor_family = A2_noncyclic_boundary_sign_stabilizers`
- `inherits_from_ce6_v1_loop = false`
- `legacy_bridge_used = false`

## Obowiązki dla PC2

1. Jawny map-anchor bez cichego selector-closure.
2. Z góry zamrożone metryki (`sign_flip_rate`, `selector_drift`) i progi.
3. Nowy baseline artifacts (oddzielne od klasy `P1388/P1389/P1390`).
4. Obowiązkowy robust-pass na perturbacjach przed jakimkolwiek exportem `L-B1-03`.

## Kluczowe wyjaśnienie: QW-2191 w R8.1 vs teraz

### Co domknęło się w Release 8.1 (QW-2191)

- Domknięto **warstwę operacyjną strict**:
  - stabilizacja pipeline’u,
  - rozdział legacy/strict,
  - governance (`no silent transfer`, rollback, audit gates),
  - gotowość do uczciwego PASS/FAIL na artefaktach.
- To jest domknięcie metodologiczno-operacyjne, nie pełny theorem-level bridge do `L_SM + L_GR`.

### Co próbujemy domknąć teraz

- Domknięcie **theorem-level** konkretnego brakującego ogniwa B1 (`L-B1-03` / `ce6`).
- To wymaga globalnej odporności i exportu dowodu (albo formalnej przeszkody),
  a nie tylko działającego pipeline’u.
- `P1391` pokazał, że obecna klasa providerów nie domyka tego warunku (obstruction v1),
  więc teraz uczciwie zmieniamy klasę providerów (`PC2`) zamiast powtarzać pętlę.

## Werdykt kroku

`PC2_STATUS := BASELINE_FROZEN_READY_FOR_FIRST_RUN`

## Następny krok

`P1393_STRICT_ONLY_PC2_FIRST_SIGN_STABILITY_AND_SELECTOR_DRIFT_RUN`

## Omówienie dla laika

R8.1 domknęło „jak testować uczciwie”.
Teraz próbujemy domknąć „czy matematycznie naprawdę działa”.
To dwie różne rzeczy: pierwsza to dobra procedura, druga to trudny dowód.
