# P1689 S639 Strict Spin-2 Curved Background Operator Sign Witness Packet (PL)

Status: `P1689_EXECUTED_STRICT_SPIN2_OPERATOR_SIGN_WITNESS_NO_FALSE_PASS`  
As of: `2026-05-14`

## Cel

Następny uczciwy krok po P1688: przejść z poziomu proxy na poziom
operatorowego testu znaku dla sektora spin-2 na tle zakrzywionym, bez odwołań
legacy i bez deklaracji fałszywego domknięcia.

## Kontekst łańcucha

`K_strict -> współczynniki -> pełny L_total -> EOM -> operator fluktuacji spin-2`

z ograniczeniem strict-core: brak mostu legacy i brak automatycznego PASS dla
renormalizacji/unitarności/background-independence.

## Co robi P1689

1. buduje minimalny operator spin-2 wokół stałokrzywiznowego tła,
2. wyznacza wartości własne operatora dla zadanej siatki modów,
3. wykonuje test znaków widma (`min_eigenvalue > 0`) jako warunek konieczny
   lokalnej stabilności/unitarności kwadratowej,
4. raportuje wynik jako `KEEP_OPEN` jeśli warunek nie domyka całej QG.

## Wynik

- test znaków widma jest przejrzysty i replikowalny,
- brak fałszywego PASS: nawet przy dodatnim widmie lokalnym pozostają otwarte
  obowiązki globalne (kontrterminy, BRST/unitary cuts, background independence
  kwantowe).

## Release context

Repozytorium jest po `release 8.1`; ten krok dopina stricte techniczny tor
`8.2/8.3` jako artefakt DOI-ready, ale bez marketingowej deklaracji domknięcia
ToE.

## Dla laika

To jak test zawieszenia w aucie: sprawdzamy, czy przy małych drganiach układ
sam się nie „rozbiega”. Jeśli test lokalnie jest dobry, to super — ale to nie
znaczy jeszcze, że auto jest gotowe na każdy teren. Nadal trzeba sprawdzić
trudniejsze warunki kwantowe.

## Następny uczciwy krok

Zbudować strict-core eksport kontrterminów 1-loop dla sektora spin-2/SM mix na
tym samym tle i połączyć go z testem unitarności cięć (Cutkosky/BRST) w jednym
witnessie PASS/FAIL/KEEP_OPEN.
