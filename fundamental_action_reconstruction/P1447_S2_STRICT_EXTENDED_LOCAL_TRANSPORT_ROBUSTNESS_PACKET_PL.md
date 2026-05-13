# P1447 S2 Strict Extended-Local Transport Robustness Packet (PL)

Status: `P1447_S2_PROTOCOL_READY_STRICT_ONLY_NO_LEGACY_BRIDGE`
As of: `2026-05-13`

## Cel

Wykonać kolejny uczciwy krok po S1/P1446:

```text
S2 = extended-local transport robustness
```

czyli sprawdzić, czy lokalny witness utrzymuje się przy kontrolowanym transporcie parametrów w klasie `kernel-split-robust`.

## Kontrakt S2

Wejście: zestaw scenariuszy transportowych `scenario_i` z polami:
`delta_scale`, `delta_phase`, `margin_before`, `margin_after`, `replay_gap`.

Warunki pass (lokalne, nieglobalne):

1. `margin_after >= margin_before` dla każdego scenariusza,
2. `replay_gap <= replay_tol` dla każdego scenariusza,
3. `strict_core_qw2191_closed == false` i `scope_of_pass == LOCAL_ONLY_NON_GLOBAL_CLAIM`.

## Semantyka werdyktu

- `PASS_LOCAL_TRANSPORT_ONLY`
- `FAIL_TRANSPORT_MARGIN`
- `FAIL_TRANSPORT_REPLAY`
- `FAIL_SCOPE`

## Decyzja profesorska

S2 jest właściwym krokiem, bo nie udaje zamknięcia globalnego, a jednocześnie zwiększa wymagania względem S1: przechodzimy z punktowego testu monotoniczności do rodziny perturbacji transportowych.

## Rekomendacja następnego uczciwego kroku

**Uruchomić P1447 na prerejestrowanym zestawie scenariuszy i dopiero po PASS przejść do S3 (minimalnego certyfikatu kompatybilności łańcucha S1+S2, nadal local-only).**

## Omówienie dla laika

W S1 sprawdziliśmy „czy działa w jednym miejscu”. W S2 sprawdzamy „czy nadal działa, gdy lekko przekręcimy pokrętła”.
Jeśli system psuje się od małych zmian, to nie jest jeszcze gotowy na większe wnioski.
