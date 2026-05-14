# P1668 / S618 Kernel <-> EOM Bidirectional Constraint-Match Packet

Status: `P1668_EXECUTED_KERNEL_EOM_BIDIRECTIONAL_MATCH_NO_FALSE_PASS`
As of: `2026-05-14`

## Cel

Wykonać następny uczciwy krok strict-only: połączyć diagonalizację proxy (`S617`) z
jawnie zapisanymi constraintami EOM i sprawdzić dwukierunkową zgodność:

`kernel -> EOM constraints` oraz `EOM constraints -> kernel relations`.

## Zakres

- forward: z parametrów kernela budujemy współczynniki `(Mpl2,cR2,cRic2,cRiem2,muH2,lambdaH,xiHR)`,
- forward: z nich generujemy proxy-constraints EOM,
- reverse: z constraintów odtwarzamy `(beta,A,omega,eta)` i relacje GR,
- mierzymy błędy rekonstrukcji oraz residuale constraintów.

## Polityka statusu

Lokalny PASS możliwy tylko dla rekonstrukcji numerycznej.
Globalny status pozostaje `OPEN_OBLIGATION` do czasu theorem-level dowodów
renormalizacji, unitarności i background-independence.

## Rekomendowany następny uczciwy krok

`S619`: zastąpić proxy-constraints jawnie wyprowadzonymi równaniami ruchu z wariacji pełnego `L_total` (metryka + pola SM) i dołączyć testy na wielu punktach parametrycznych strict-kernel.

## Omówienie dla laika

To etap „sprawdź obie strony równania”: czy z ustawień startowych dostajesz
te same zasady ruchu, i czy z tych zasad da się wrócić do ustawień.
To buduje zaufanie do spójności modelu, choć nie zamyka jeszcze całej ToE.
