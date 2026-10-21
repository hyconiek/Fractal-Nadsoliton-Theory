# Długi czas pamięci: opłacone oszacowanie oscylacyjne, nie prawo tłumienia

Data: 2026-10-21. Nowy warunkowy lemat operatorowy dla rzeczywistego
jądra NM. Bez nowego prawa, szumu, danych wejściowych lub niezależnej recenzji.

## Wynik

Udowodniłem przybliżenie przyczynowej pamięci dla zadanych wolnych obwiedni
aż do czasu S epsilon^-2. Mimo rozbieżnego pierwszego momentu bezwzględnego
kontinuum, w faktycznym nierezonansowym paśmie źródła błąd wynosi
O(epsilon^2), także po jednym różniczkowaniu w oryginalnym czasie.

Nie zastępuje się jednak całej pamięci samą podatnością. Jawnie pozostaje
wyraz początkowy, którego część pochodząca od modu kształtu nie zanika.
Dla stałej obwiedni jest to dokładne rozwiązanie przyczynowe, nie dowolny
stan ustalony wstawiony w miejsce zerowych danych Cauchy'ego.

Metoda korzysta z ogona oscylacyjnego jądra, nie z fałszywego założenia
integral t|k_c(t)|dt<infinity. Ten ogon jest całkowalny poza częstotliwościami
progowymi; dla modu kształtu użyłem jawnych ograniczonych pierwotnych.
Wolna siła o amplitudzie epsilon^2 ma w tej projekcji błąd O(epsilon^4).

## Zachowane wyniki przeciwne

Na rzeczywistym progu kontinuum przybliżenie traci rząd: dla obwiedni
liniowej błąd jest rzędu epsilon, nie epsilon^2. Ponadto jednostajna
kontrola tylko pierwszej pochodnej obwiedni nie wystarcza nawet przy
dozwolonej bazowej częstości0. Mała szybka składowa może pobudzić mod
kształtu rezonansowo i dać błąd rzędu jedności. Te kontrprzykłady pokazują,
dlaczego nie wolno pomijać luk częstotliwości i drugiej pochodnej.

## Co faktycznie zamknięto

Opłacono część operatorową kroku2 w NM/NEXT_OPEN_LEMMA, z zachowaniem
koherentnej pamięci początkowej. NIE udowodniono jeszcze, że rzeczywista
amplituda pełnego pola ma wymaganą jednostajnie C2 wolną reprezentację.
NIE przedłużono błędu pełnego PDE na epsilon^-2 ani nie udowodniono
uśrednionego równania zaniku. Rozróżnienie jest zapisane w
[GATE_DELTA.md](GATE_DELTA.md).

To potrzebny nowy składnik wyprowadzenia efektywnej dynamiki z obecnych
sprzężonych równań, a nie nowe równanie mikroskopowe lub źródło losowości.
Nie ma wniosku o fizycznym czasie, energii, obserwatorze, QM/GR lub ToE.
F01 pozostaje otwarte. Stare źródła, plany i AGENTS zachowano.

[Pełny dowód i kontrprzykłady](PROOF_OR_OPEN_LEMMA.md),
[kontrole rachunkowe](RESULTS.json), [odtworzenie](REPLAY.md).
Przechodzi 41 nowych ograniczonych kontroli OM oraz osobno odtworzenie
34 wcześniejszych kontroli NM. Związano 12 wejść i zachowano 35 wpisów
manifestów. Nie są to niezależne dowody ani symulacja PDE na długim czasie.
