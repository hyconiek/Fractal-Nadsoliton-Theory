# P1997 / S947 — Strict Cutkosky channelwise state-sum solver witness

P1997 wykonuje następny uczciwy krok po P1996: zamiast globalnego profilu korekty
cięcia eksportuje jawny, skończony sumator kanał-po-kanale dla stanów pośrednich
(`gg`, `gh`, `hh`) i liczy `Delta_opt(s)=ImM-CutSum` na tej samej siatce energii.

## Zakres

- strict lane only,
- model proxy dla skończonej klasy stanów pośrednich,
- brak claimu theorem-grade closure pełnego UR-link.

## Wynik

- Eksportowano jawne termy kanałowe `Cut_c(s)` oraz ich udział normowy.
- `Delta_opt` pozostaje niezerowe (kontrolowane), więc obiekt rozróżnia: 
  *niedomknięcie sumy stanów* vs *potencjalna obstrukcja strukturalna*.
- Zachowano guardrail `OPEN_OBSTRUCTION_WITH_TRACE` i brak fałszywego domknięcia.

## Następny krok

Kalibracja amplitud kanałowych z jawnych wkładów pętlowych i rozszerzenie listy
stanów pośrednich tak, by test `Delta_opt` był oparty o rzeczywistą sumę cięć,
a nie proxy-parametryzację.
