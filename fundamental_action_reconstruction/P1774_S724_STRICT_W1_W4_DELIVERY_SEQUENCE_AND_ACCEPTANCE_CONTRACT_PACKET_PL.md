# P1774 / S724 — Sekwencja dostaw W1–W4 i kontrakt odbioru

Status: `P1774_S724_W1_W4_ACCEPTANCE_CONTRACT_READY_NO_FALSE_PASS`
As of: `2026-05-15`

## Technical progress

Po synchronizacji stanu (`P1773`) wprowadzono twardy kontrakt odbioru dla
bloków W1–W4, żeby przejść z planu do wykonania bez rozmycia kryteriów jakości.

Zdefiniowano:

1. sekwencję dostaw: `W1 -> W2 -> W3 -> W4`,
2. kryteria odbioru każdego bloku jako `FULL_EXPORT`,
3. warunki globalne przed ponownym uruchomieniem `G_BW`.

## Co zostało dowiedzione

1. Samo „dostarczono wzór” nie wystarcza — każdy Wi musi mieć:
   - jawny komponentowy symbol list,
   - lock konwencji indeksowej,
   - mapę projekcji do `B1/B2/B3/C1/C2`.
2. To eliminuje ryzyko pozornego postępu typu scaffold-only.

## Co nadal jest OPEN

1. Rzeczywiste dowiezienie W1..W4.
2. Kolejne uruchomienie `G_BW` po odbiorze wszystkich Wi.
3. Odblokowanie BRST/Cutkosky.

## Ryzyka false-pass

1. Oznaczenie Wi jako „gotowe” bez klasyfikacji `FULL_EXPORT`.
2. Zmiana tła lub bazy residualu między iteracjami.
3. Uruchomienie BRST/Cutkosky przed końcowym werdyktem `G_BW`.

## Następny uczciwy krok

Dostarczyć `W1_H_R2_componentwise` jako pierwszy pełny eksport wraz z mapą
projekcji do `B1/B2/B3/C1/C2` i dopiero wtedy przejść do `W2`.

## Dla laika

To jak lista odbiorowa w inżynierii: każdy moduł musi przejść ten sam zestaw
kontroli jakości, zanim pozwolimy systemowi iść dalej.
