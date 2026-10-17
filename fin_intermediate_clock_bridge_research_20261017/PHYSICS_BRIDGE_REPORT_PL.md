# DH po audycie: most do fizyki wymaga także zegara i instrumentu

Data: 2026-10-17. Analiza teoretyczna, nie pomiary, formalizacja ani
niezależna recenzja naukowa. Poprawne wyniki samego DH są przyjęte osobno
w [intake](../fin_dense_orbit_fibre_transport_review_20261017/INTAKE_20261017.md).

## Co nowego daje DH

Udowodniono Hall dla jawnej rodziny gęstych przekrojów: usuwa się dwanaście
dobrze rozdzielonych stanów jednej orbity. Dodatnie marginesy prądu usuniętej
orbity dominują masę ewentualnych brakujących kolumn na prawdziwym brzegu.
To rzeczywisty postęp poza rzadkimi zbiorami z poprzedniego PC. Nie jest
jednak twierdzeniem dla dowolnych skupionych dopełnień ani pełnym sprzężeniem.
Nawet pełny Hall nie dostarczyłby sam źródła fizycznego nośnika lub kinetyki.

## Rzeczywiście wykonana lektura

W działającym lokalnym browser-harness wyszukałem prace o eliminacji
intermediatów, odczytałem metadane Cappellettiego–Wiufa i pełną stronę
nowej pracy Hoessly'ego–Wiufa–Xia: online2025, wydanie2026. Odpowiednie
sekcje PDF-ów oraz erratum Esposita czytano dodatkowym narzędziem web,
gdy Firecrawl zwrócił HTTP402. Zakresy i ograniczenia są w
[bazie literatury](LITERATURE_DATABASE.md). Nagrywanie pozostało wyłączone;
bez kupowania kredytów, usług w chmurze lub sprzętu.

Literatura odróżnia redukcję stanów, zbieżność w ustalonych chwilach,
pełne prawa trajektorii i stacjonarność. Żaden z tych wyników nie przenosi
automatycznie pozostałych. Ujawnione erratum przypomina też, że warunki
termodynamiczne dotyczące liczby rezerwuarów trzeba czytać dosłownie.

## Nowy most warunkowy: wspólny pośrednik zamiast wpisanego normalizatora

Wyprowadziłem [porównawczy model z jedną aktywną pustą kopią](TRACE_CLOCK_AND_INSTRUMENT_THEOREM.md).
Widoczny stan aktywuje wybraną kopię z szybkością nu; pośrednik uzupełnia
etykietę b z szybkością kappa exp[(g/N)(Ar)_b]. Z praw konkurujących
wyjść powstaje dokładnie softmax. Po wymazaniu czasu w pośredniku proces
śladowy jest **dokładnym pełnym native FIN**, z poprawnymi licznościami.
Rozszerzony proces jest odwracalny i ma native Gibbs po warunkowaniu na
stany widoczne. To zastosowanie znanej konstrukcji śladu, nie odkrycie
ogólnej procedury Schura i nie wyprowadzenie fizycznej pustej kopii.

Ten proces nadal ma natywne przejścia F. Nie jest poszukiwanym sprzężeniem
diagnostycznym na H i nie zamyka DH-002.

## Gdzie pojawia się fizyka, a gdzie tylko zmiana opisu

W czasie rzeczywistym jeden cykl trwa Exp(N nu)+Exp(kappa Z), a nie samo
Exp(N nu). Wymazywanie czasu pośrednika nie jest kalibracją zegara.
Dla jawnie zdefiniowanego odczytu „ostatni kompletny stan” zapłaciłem
ograniczenie pełnego prawa ścieżki:

    TV <= min{1,(N nu T+1) N nu/(N nu+kappa c_Z)},
    c_Z=12 exp(-|g| max|A_ij|)>0.

Na horyzoncie krytycznym sqrt(N)S/nu wystarcza kappa_N/nu >> N^(5/2).
To konserwatywny warunek wystarczający **dla tej szeregowej konstrukcji**,
nie konieczna granica, optimum, praktyczny budżet ani stała natury.

Dla surowego odczytu wykrywającego pośrednik TV jest natomiast co najmniej
1-exp(-N nu T), niezależnie od kappa. Krótki czas zajęcia nie oznacza
niewidoczności w idealnym ciągłym zapisie. Fizyczny instrument musiałby
osobno uzasadnić cenzurowanie/rozdzielczość; nie wolno usuwać „invalidów”
po obejrzeniu danych. Także odbicie w D_N trzeba wykonywać we właściwej
kolejności: odrzucanie wyjść przed eliminacją pośrednika zmienia normalizator.

## Konkretny brakujący obiekt fizyczny

[Analiza barier](ACTIVATION_BARRIER_SOURCE_GATE.md) pokazuje, czego wymaga
mechanizm resetu: jednakowej **względnej** bariery aktywacji B_a-E_a,
odpowiednich prefaktorów, licencjonowanego pośrednika i jego szybkiego zegara.
Wspólna bariera **absolutna** też daje poprawny Gibbs, ale inną kinetykę.
Samo równowagowe oddziaływanie A7 nie wybiera między tymi mechanizmami.
Wyliczenie bariery z działania opartego już na HB byłoby kołem logicznym.

Nie znaleziono niezależnego źródła tych barier z obecnych równań FIN.
Porównawczy model pokazuje **jak mógłby wyglądać warunkowy most**, nie
dlaczego natura miałaby go realizować. AS/KB/QP pozostają nieopłacone.

| Warstwa | Nadal potrzebny obiekt |
|---|---|
| Nośnik i A7 | Niezależna od celu selekcja stanów, wsparcia i wag |
| Kinetyka | Pochodzenie barier/prefaktorów, nie samo detailed balance |
| Czas i energia | Fizyczne nu,kappa,J,beta i bilans całego układu |
| Kompozycja/lokalność | Źródło incydencji i propagacji; kopie nie są przestrzenią |
| Odczyt | Uzasadniony instrument, rozdzielczość, zakłócenie i invalidy |
| Empiria | Niezależne mapowanie, kalibracja, custody i walidacja |
| QM/GR/ToE | Własne operacje nieklasyczne/geometria i fizyczne mosty, nadal brak |

## Następne badania

W matematyce warto przejść do **jednego zamrożonego klastra o dodatniej
masie granicznej**, nie mnożyć skończonych rozdzielonych wzorów. Intake
wyprowadza geometrię erozji takiej kuli; cienka warstwa siatki i kontrola
różnicy sąsiednich włókien pozostają zadaniami.

W fizyce wymagana jest jedna rzeczywiście nowa przesłanka/ewolucja,
która niezależnie źródłuje bariery i nośnik. Jeśli jej nie ma, zatrzymać
promocję fizyczną, nie przerabiać kolejnego wygodnego liftu na źródło świata.
[Następny łańcuch](../FIN_MESOSCOPIC_HALL_AND_INTERMEDIATE_SOURCE_NEXT_CHAIN_20261017.md)
rozdziela te dwie bramki. Nie dodano nowego prawa do primary FIN ani
informacyjnej warstwy pod nadsolitonem. Wynik nie jest QM/GR/SM/ToE.
