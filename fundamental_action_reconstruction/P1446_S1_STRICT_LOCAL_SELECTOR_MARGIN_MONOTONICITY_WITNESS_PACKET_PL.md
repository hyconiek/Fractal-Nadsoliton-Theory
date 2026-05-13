# P1446 S1 Strict-Local Selector Margin Monotonicity Witness Packet (PL)

Status: `P1446_S1_PROTOCOL_READY_STRICT_ONLY_NO_LEGACY_BRIDGE`
As of: `2026-05-13`

## Cel

Wykonać pierwszy twardy krok po P1445:

```text
S1 = strict-local selector margin monotonicity witness
```

na trasie `kernel-split-robust`, bez globalnych claimów i bez transferu ról legacy.

## Kontrakt testu S1

Dla ustalonej siatki perturbacji `eps_grid` oraz wektora marginesów selektora `margin(eps)`:

1. monotoniczność: `margin(eps_{i+1}) >= margin(eps_i)` dla każdego `i`,
2. replay A/B: `max |margin_A(eps)-margin_B(eps)| <= replay_tol`,
3. pass-scope: wynik `PASS_LOCAL_ONLY` nie może być opisany jako `GLOBAL_CLOSURE`.

## Semantyka wyników

- `PASS_LOCAL_ONLY`: spełnione 1-3, ale bez global closure claim.
- `FAIL_MONOTONICITY`: złamany warunek 1.
- `FAIL_REPLAY`: złamany warunek 2.
- `FAIL_SCOPE`: złamany warunek 3.

## Decyzja profesorska

Najbardziej uczciwa decyzja metodologiczna to uruchomienie S1 na minimalnym, jawnie opisanym protokole z eksportem JSON i bez narracyjnej inflacji znaczenia PASS.

## Rekomendacja następnego uczciwego kroku

**Uruchomić checkpoint P1446 i opublikować artifact `p1446_s1_selector_margin_monotonicity_summary.json` jako jedyne źródło werdyktu S1.**

## Omówienie dla laika

To jest jak test stabilności mostu w małej skali: jeśli mały odcinek nie jest stabilny i powtarzalny, nie wolno mówić, że cały most jest gotowy.
Jeśli mały test przejdzie, to nadal nie znaczy „wszystko gotowe” — tylko że warto iść do kolejnego, trochę większego testu.
