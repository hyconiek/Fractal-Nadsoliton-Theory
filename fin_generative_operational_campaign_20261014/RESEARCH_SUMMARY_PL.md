# FIN — wyniki kontynuacji RG-E01–RG-E05 (14.10.2026)

## Zakres i ostateczny status

Kampania wykonała lokalne dowody i obliczenia skończonego modelu, ale **nie przeprowadziła eksperymentu fizycznego**. Zgodnie z instrukcją RG-E06 pozostaje do decyzji architekta. Status: **DESIGN_ONLY / NO_NEW_SOURCE / NO_NEW_PHYSICAL_INPUT**. Wszystkie liczby niżej dotyczą z góry zadanych matematycznych modeli.

## RG-E01: warstwa operacyjna

Rozdzielono M3: trzy etykietowane kopie q12, g=5.145228719489142, obserwacja dwóch pierwszych, od projektu PHYS-014: N2, dwa węzły q12, g=3, obserwacja wyłącznie różnicy etykiet. PHYS-014 jest projektem pomiarowym; nie udostępniono rzeczywistych danych kanałów, zegara, błędów detektora czy kalibracji. N2 nie waliduje pamięci M3.

## RG-E02: identyfikowalność i swobody cechowania

Dla rzeczywistych rozkładów warunkowych M3 wyprowadzono dokładną tożsamość logarytmicznych ilorazów:

D(j,l;b,b';c) = (g/3) [A_jb-A_lb-A_jb'+A_lb'].

Jest ona niezależna od zegara, wspólnej etykiety c i stałego pola jednociała. Pełne kontrasty wyznaczają tylko podwójnie scentrowane gA — nie oddzielają g od skali A; nie rozstrzygają także o fizycznym pochodzeniu macierzy. Sam projekt N2 nie daje dostępu do takiego warunkowania M3.

**Teoretyczne współczynniki Fouriera FIN** na aktywnych modach k=3,4,5,6, w proporcji do k3: 1,000000; 1,121424; 1,171917; 1,194134. Płaski P7 ma wszystkie równe 1. To test porównania zadanych kształtów, nie dowód konieczności rangi siedem.

## RG-E03: separacja modeli N2 przy g=3

| Porównanie idealnych rozkładów różnicy etykiet | TV |
|---|---:|
| FIN / full Potts trace | 0,0969987870350775 |
| FIN / flat P7 trace | 0,0117871915793791 |
| Full Potts / flat P7 | 0,1063401407356514 |

Dla tych trzech porównań przeprowadzono *osobne* obliczenia przedziałowe exact-decimal (założona poprawność `mpmath.iv`) i otrzymano przedziały o szerokości ~1e-15. Te wyniki nie posiadają automatycznie pełnego certyfikatu błędów fizycznego toru odczytowego.

Przy hipotetycznych niepewnościach ±1% parametru g i TV detektora 0,005 na model, konserwatywny wzór TV≥TV_nominal-r1-r2 **nie daje dodatniej gwarancji** oddzielenia FIN od flat P7. Nie jest to dowód rzeczywistego zerowego dystansu przy tych konkretnych błędach. Bez jakiegokolwiek ograniczenia detektora istnieje dokładny kontrprzykład (kanał całkowicie zacierający rozróżnienie), więc minimalna separacja wynosi zero.

## RG-E04: kinetyka, nie tylko równowaga

Dla osobnego N2 porównano heat-bath, Metropolisa i Barkera. Po dobraniu *modelowych* częstości prób tak, aby przy g=0 każdy generator dawał taką samą relaksację pojedynczej kopii, we wszystkich trzech przypadkach zachowano ten sam rozkład Gibbsa przy g=3.

Dla jednolicie przygotowanych dwóch węzłów (g=0), skoku sprzężenia do g=3 oraz idealnego dwuczasowego odczytu różnicy przy t=0,4 i t=1,0 wyliczono pełne tablice 13×13, z dodatkowym wynikiem INVALID. W modelu FIN odległości TV wynoszą:

- heat-bath kontra Metropolis: **0,2566279210**;
- heat-bath kontra Barker: **0,1948765862**;
- Metropolis kontra Barker: **0,0850564178**.

Są to różnice teoretyczne w czasie generatora. Bez niezależnego pomiaru relaksacji g=0 na aparaturze nie uzasadniają identyfikacji jej prawdziwego zegara, a bez charakterystyki detektora nie gwarantują separacji eksperymentalnej.

## RG-E05: zamrożony cel projektu / stop-certyfikat

**Test pierwotny:** dwuwęzłowy N2 g=3, 12-binowy histogram różnicy etykiet, FIN przeciwko flat P7 (przy skontrolowanych pozostałych modelach jako alternatywach).

**Test wtórny tylko po kalibracji zegara:** 13×13 dwuczasowy rozkład różnicy etykiet, heat-bath kontra Metropolis, przy tych samych regułach przygotowania i odczytu.

Przed rzeczywistą walidacją wymagane są niezależnie ustalone i zapisane: identyfikatory 24 kanałów, wagi i temperatury, kalibracja g=0, błędy/invalid detektora, stabilność, zależność w czasie i efektywna liczba niezależnych obserwacji, rozdzielenie kalibracji i walidacji. Nie podjęto zakupu ani kontaktu z laboratorium, nie pobierano nowych danych.

## Sprawdzalność pakietu

W zestawie są pełne wektory 12-binowe, macierze 144×144 implicit/sparse i 13×13 wynikowych praw dwuczasowych, trzy skrypty badawcze, niezależny weryfikator, kod raportowania, manifest SHA-256 oraz kompletny poprzedni handoff RG-D. **58/58 testów niezależnego weryfikatora** i pięć identycznych bajtowo wybranych wyjść przy ponownym wykonaniu w osobnym katalogu. Obliczenia dynamiki pozostają numeryczne, bez formalnego pełnego oszacowania błędu wszystkich operacji zmiennoprzecinkowych.

**Nie udowodniono:** niezależnego źródła A7, fizycznej realizacji konkretnego generatora, lokalności, pola, cząstek, fizycznego czasu, mechaniki kwantowej, grawitacji ani teorii wszystkiego.
