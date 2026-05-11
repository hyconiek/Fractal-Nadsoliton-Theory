# P1147 Strict Selector Obstruction Location Probe Packet

Status: `P1147_EXECUTED_STRICT_SELECTOR_OBSTRUCTION_LOCATION_PROBE_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Wykonać bardziej ścisły następny krok po `P1146`: nie tylko stwierdzić
oscylację, ale zlokalizować **pierwszą przeszkodę znaku** i sprawdzić, czy da
się ją usunąć przez dodatnie Shannon-like ważenie.

## Candidate family

Badana rodzina strict-side:

```text
S_gamma(d) = exp(-alpha*d - gamma*d^2) * cos(omega*d+phi)/(1+beta*d^eta), gamma>=0
```

z parametrami strict tuple:

```text
alpha=4 ln 2, omega=0.18575, phi=0.16250, beta=1, eta=1.8
```

## Professor-level decision

Wybieram test „obstruction-location”:

1. analitycznie wyznaczyć kandydat na pierwszy zero-crossing z fazy
   `cos(omega*d+phi)`,
2. numerycznie potwierdzić pozycję pierwszej zmiany znaku dla kilku `gamma`,
3. sprawdzić niezmienniczość tej przeszkody względem dodatnich wag.

To jest uczciwy krok, bo testuje strukturę przeszkody, a nie tylko pojedynczy
przypadek siatki.

## Result

Wynik z artefaktu maszynowego:

1. pierwsza przeszkoda znaku jest praktycznie stała dla badanych `gamma`,
2. lokalizacja zgadza się z analityczną predykcją fazową,
3. dodatnie monotoniczne ważenie nie usuwa pierwszego flipu znaku,
4. audit proxy pozostaje `QW-2191: BLOCKED`.

## Artifact

- probe script:
  `p1147_strict_selector_obstruction_location_probe.py`
- generated summary:
  `generated/p1147_strict_selector_obstruction_location_probe_summary.json`

## Honest boundary

`P1147` wzmacnia rygor przez wskazanie, że blokada nie jest artefaktem doboru
jednej siatki ani jednej wagi, tylko wynika z fazowej struktury kandydata.

Nie jest to jeszcze theorem-level no-go dla wszystkich możliwych strict
konstrukcji selektora.
