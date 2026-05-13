# P1391 Strict-only `ce6` Formal Obstruction Theorem Export (No Legacy Bridge) — PL

Status: `P1391_EXECUTED_FORMAL_OBSTRUCTION_EXPORT_V1`
As of: `2026-05-13`

## Cel

Po `P1390` (robust-fail) formalnie wyeksportować theorem-level obstruction dla `ce6`,
aby zamknąć obecny cykl lokalnych prób i utrzymać dyscyplinę no-false-pass
w torze strict-only.

## Rygor

- `legacy_bridge_used = false`
- brak transferu ról legacy->strict
- jawna kwalifikacja: obstruction theorem, nie closure theorem

## Treść eksportu (v1)

`OBSTR-CE6-v1`:
Dla obecnej klasy providerów oraz boundary-weight map `W_boundary^(v1)`
nie istnieje stabilny uniform bound gwarantujący
`sign_flip_rate <= epsilon_sign_v1` na całej rodzinie dopuszczalnych perturbacji.

## Podstawa dowodowa

- `P1388`: poprawa lokalna, ale threshold nieosiągnięty.
- `P1389`: adaptacyjna re-waga poprawia wynik, ale nadal nad progiem.
- `P1390`: sweep robustności daje `ROBUST_FAIL` (`pass_count < total`).

## Werdykt

`CE6_FORMAL_STATUS := OBSTRUCTION_THEOREM_EXPORTED_V1`

Konsekwencje:
- `L_B1_03_EXPORT_STATUS = NOT_EXPORTED`
- `B1 = OPEN`
- zakończenie bieżącej klasy prób lokalnych

## Decyzja profesorska

Następny krok: `P1392_STRICT_ONLY_NEW_PROVIDER_CLASS_ANCHOR_DESIGN`
- zaprojektować nową klasę provider/anchor (noncyclic),
- uruchomić nową gałąź z odseparowanym baseline,
- nie wracać do powtarzania tych samych loopów patch-weight w tej klasie.

## Omówienie dla laika

To jak oficjalny raport: obecna technika poprawy nie daje stabilnego wyniku.
Nie brniemy w nieskończone poprawki tego samego typu,
tylko uczciwie zamykamy etap i przechodzimy do nowej metody.
