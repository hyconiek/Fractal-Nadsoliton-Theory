# Nowe badanie F01: brakującego prawa nie zastępuje jeden średni potencjał

Data: 2026-10-20. Badanie wykonane samodzielnie, bez subagentów, zmiany
AGENTS lub oryginalnych pakietów. Dowody warunkowe względem przyjętych
matematycznych danych FIN, nie formalizacja ani niezależna recenzja.

## Co rzeczywiście wyprowadziłem

**Dokładne wymagane pole historii.** Dla pełnych zadanych wag pi i legalnego
nośnika K istnieje jedno pole f=log P(x,x), z dokładnością do stałej, które
wybiera te wagi przez maksymalizację entropii historii. Jest to obiekt
matematyczny obliczony z celu, nie niezależne źródło tego celu.

**Konkretną przeszkodę dla uproszczonego mostu.** Dla rzeczywistego A7 i
dwóch oznaczonych kopii pole to ma dodatnie harmoniczne 1 i 2, choć sam
potencjał oddziaływania ich nie ma. Dowód obejmuje każde g>0. Zatem żadna
zmiana jednego mnożnika przy średnim oddziaływaniu nie odzyska pełnej
zadanej miary Gibbsa. Nie trzeba opierać tego wniosku na numerycznym
poszukiwaniu minimum ani uprzednio zakładać Markowowskości historii.

**Jawnego konkurenta.** Zmiana p_j na p_j-t cos(2pi j/12) zachowuje
normalizację, średnie oddziaływanie i symetrie, ale zwiększa dostępną
entropię historii. Przy zamrożonym g=12/Lambda6 opłaciłem dodatniość i
krzywiznę na całym potrzebnym przedziale. Dla dokładnego t=10^-16 wzrost
jest większy niż 1/(2*10^18) nata na oznaczony krok. To konserwatywny
świadek istnienia, nie praktyczna lub zmierzona wielkość fizyczna.

[Dowód pola](INVERSE_FIELD_THEOREM.md), [rzeczywisty A7](A7_MISSING_HARMONIC_THEOREM.md),
[konkurent i oszacowanie skończone](EXPLICIT_SAME_MEAN_ALTERNATIVE.md).

## Dlaczego to ma znaczenie dla źródła ewolucji

Poprzednia konstrukcja skończonych historii wymagała całej miary Gibbsa.
Sprawdziłem teraz konkretną możliwość usunięcia tego założenia. Wynik ją
wyklucza w zamrożonej klasie: **jeden warunek średniego oddziaływania to
za mało**. Potrzebne jest dodatkowe niezależne ograniczenie historii lub
wag. Samo policzenie wymaganego pola f z założonej pi byłoby rozumowaniem
kołowym jako źródło F01.

Nie jest to kontrprzykład do FIN-HYP-HIST-ME-001 z CAŁĄ miarą pi, do HB
ani do każdej możliwej teorii źródłowej. Po usunięciu oznaczonych
multiplikacji i zmianie nośnika na kompletny graf separacji wynik znika:
tam pole log pi wystarcza. Ta kontrkontrola pokazuje istotność reprezentacji;
nie podmieniam nią zamrożonego modelu.

## Sprawdzenie i decyzja

Weryfikator ma **42 ograniczone kontrole**: rachunek wymierny i Q(sqrt(3)),
normalizacje Fouriera, dokładne przedziały spektralne A7, szereg z ogonem,
certyfikat stałych oraz osobny 16-stanowy przykład pochodnych entropii.
Wiąże 16 wejść i zachowuje 10 manifestów/161 wpisów poprzedników. Te kontrole
nie są 42 dowodami asymptotycznymi. N2 przy wartości g_c nie udaje obszaru MT.
[Certyfikat](CERTIFICATE.json), [wyniki](RESULTS.json), [odtworzenie](REPLAY.md).

Wagi Gibbsa, wybór nośnika i zasada wyboru historii nadal wymagają naturalnego
uzasadnienia. **F01/F03/F06 pozostają OPEN_SOURCE.** Nie dopisuję nowego
postulatu, strojenia mnożników ani automatycznej kampanii harmonicznych.
Następny obiekt źródłowy musi wnieść realną niezależną przesłankę, a nie
nazwać wyprowadzony z celu f prawem natury. [Handoff](HANDOFF.md).
