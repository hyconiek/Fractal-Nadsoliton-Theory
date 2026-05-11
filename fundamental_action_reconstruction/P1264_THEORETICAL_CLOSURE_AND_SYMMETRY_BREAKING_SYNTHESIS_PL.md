# P1264 — Synteza: teoretyczne domknięcie i źródło łamania symetrii

## 1) Co znaczy „teoretyczne domknięcie” w tej architekturze?

W tym repo **domknięcie teoretyczne** nie oznacza „mamy dużo pozytywnych testów”, tylko:

1. mamy formalny wynik `B1` (**bridge theorem**) **albo** `NB1` (**non-bridge theorem**),
2. mamy jawnie rozstrzygnięty problem selektora (bariera `QW-2191`) w lane strict-core,
3. nie ma ukrytego transferu ról z `K_legacy_ont` do `K_strict_gate`.

Dopiero ten zestaw pozwala przejść z „lokalna gotowość operacyjna” do
„globalna zamkniętość teorii”.

---

## 2) Co konkretnie „łamie symetrię” w tej perspektywie?

W obecnym stanie projektu **nie ma jeszcze wyeksportowanego wewnętrznego, strict-core źródła łamania symetrii**, które jednoznacznie domyka selektor.

To oznacza:
- nie wolno twierdzić, że sama dynamika strict-gate automatycznie wybiera unikalny rdzeń,
- nie wolno zamieniać „działa obliczeniowo” na „udowodnione ontologicznie”.

Praktycznie: jeśli łamanie symetrii jest użyte, musi być
- albo formalnie wyprowadzone wewnątrz strict-core,
- albo jawnie oznaczone jako dodatkowa przesłanka (nonstrict/axiom-augmented),
  bez udawania pełnej strict closure.

---

## 3) Cała teoria w skrócie (dla spójności logicznej)

### Warstwa ontologiczna
- punkt wyjścia: `nadsoliton -> light -> matter -> emergent observer`.
- historycznie-kanoniczny obiekt: `K_legacy_ont`.

### Warstwa operacyjna
- narzędzie robocze: `K_strict_gate`.
- pipeline W1 bada, czy strict-lane daje stabilne i falsyfikowalne konsekwencje.

### Warstwa pomostowa (krytyczna)
- pytanie kluczowe: czy istnieje formalny most `legacy -> strict` (`B1`),
  czy przeciwnie: formalny no-go (`NB1`).
- bez tego pytania teoria nie jest domknięta globalnie.

### Warstwa walidacyjna
- benchmarki, niepewność, replikacja (P1259-P1263) budują wiarygodność operacyjną,
  ale same nie zamykają ontologii.

---

## 4) Decyzja profesorska: co jest „następnym uczciwym krokiem”

Następny krok powinien być **teorematyczny**, nie tylko inżynieryjny:

1. sformułować minimalny kandydat `SB1` (symmetry-breaking source theorem) w strict-core,
2. wykazać jego kompatybilność z `QW-2191` (albo formalną niekompatybilność),
3. podłączyć `SB1` jako warunek do `B1/NB1` i dopiero wtedy rozstrzygać global closure.

To odcina niejawne założenia i daje uczciwą, ścisłą drogę do domknięcia teorii.
