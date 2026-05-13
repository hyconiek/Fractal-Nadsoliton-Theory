# P1501 — S4.51 QW-2191 Adversarial Falsifier Sweep (PL)

Status: `P1501_EXECUTED_QW2191_ADVERSARIAL_FALSIFIER_SWEEP_LOCAL_ONLY`
As of: `2026-05-13`

## Cel

Wykonać niezależny test odporności na obalenie dla obiektów:

- `S_internal_v1`,
- `W_Fmap_v1`.

## Decyzja profesorska

Budujemy adversarial warianty parametrów w granicach dopuszczalnych policy,
a następnie szukamy kontrprzykładu łamiącego jednocześnie:

1. dodatniość kanałów,
2. spójność orientacji,
3. redukcję luki selektora.

Brak kontrprzykładu = silne wsparcie dla global-closure candidate.
