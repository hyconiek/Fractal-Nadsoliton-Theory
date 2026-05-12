# RELEASE 8 — Podręcznik fizyki licealisty (PL)
## Domknięcie bramki O3 dla `QW-2191` w teorii nadsolitonu

**Wersja:** 8.0.1 (edukacyjna, licealna)  
**Data:** 2026-05-12  
**Zakres:** opis obecnego domknięcia strict-gate O3 w repo (`P1313..P1339`)  

---

## 0. O co chodzi w jednym zdaniu?

Chcemy pokazać, że z wielu „dopuszczalnych kompasów matematycznych” teoria wybiera **jeden spójny kierunek** i nie ma już legalnej alternatywy — to właśnie nazywamy domknięciem bramki O3 dla `QW-2191`.

---

## 1. Intuicja fizyczna (po ludzku)

Wyobraź sobie, że opisujesz Wszechświat jak układ falowy.  
Jeśli równania pozwalają na kilka równoważnych „kierunków znaku”, teoria nie jest domknięta.  
Jeśli da się uczciwie pokazać, że wszystkie legalne drogi dają ten sam wybór — wtedy bramka unikalności jest zamknięta.

W tym projekcie tę bramkę nazwano:

\[
\text{O3 Nonuniqueness Exclusion}
\]

---

## 2. Najważniejsze obiekty matematyczne

### 2.1. Kernel strict (roboczy)

W pipeline strict używany jest kernel:

\[
K_{\text{strict}}(d)=\frac{\cos(\omega d + \phi)}{1+\beta d^{\eta}}
\]

To on napędza późniejsze testy operacyjne (nie mieszamy go po cichu z legacy-kernel).

### 2.2. Wspólna reprezentacja porównawcza

Każdy kandydat selektora tłumaczymy do tej samej „tabelki”:

\[
R_{\text{common}}=(\text{branch\_sign},\;\text{phase\_class},\;\text{amplitude\_class},\;\text{residual\_slot},\;\text{admissibility})
\]

Dzięki temu porównujemy „jabłka z jabłkami”, a nie różne formaty.

### 2.3. Kluczowy score kandydata `v4`

W finalnym etapie używany był score:

\[
s(\phi,a)=(\phi-0.40)+0.9(a-0.58)-0.6(\phi-0.40)^2
\]

oraz reguła znaku:

\[
\operatorname{sign}(s)=
\begin{cases}
+1 & s\ge 0\\
-1 & s<0
\end{cases}
\]

---

## 3. Ścieżka domknięcia krok po kroku (skrót uczciwy)

1. **Katalog kandydatów** (`P1314`) — spis legalnych klas.  
2. **Normalizacja** (`P1315`) — wspólny język `R_common`.  
3. **Pairwise matrix** (`P1316`) — początkowo zostaje residual ambiguity.  
4. **Sweep + replay/adversarial** (`P1317`, `P1318`) — brak twardych kontrprzykładów, ale jeszcze nie closure.  
5. **Probe residual slotu** (`P1319`) — neutralność nie była od razu udowodniona.  
6. **Ewolucja kandydatów** (`v2`, `v3`, `v4`) — `v2` obalony, `v3` zbyt trywialny, `v4` stabilny i informatywny.  
7. **Formal obligations** (`P1327`, `P1328`) — przejście od numerics do warstwy theorem-level.  
8. **Global L2 final attempt** (`P1338`) — domknięcie ostatniego brakującego warunku.  
9. **Niezależna replikacja** (`P1339`) — wynik utrzymany na osobnym replayu.

---

## 4. Co znaczy „domknięcie” tutaj technicznie?

W checkerze O3 była lista obowiązków.  
Po `P1338` i odświeżeniu checkera wszystkie dały PASS:

\[
\text{pass\_count}=5/5,
\qquad
\text{o3\_strict\_ready}=\text{true}
\]

i status:

\[
\text{qw2191\_strict\_status}=\text{CLOSED}
\]

w sensie rygoru tego pipeline.

---

## 5. Czego to **nie** znaczy?

To NIE znaczy automatycznie:

- „wszystko w fizyce jest już domknięte na zawsze”,
- „dowolny przyszły audyt niczego nie zmieni”,
- „globalne ToE closure w najmocniejszym filozoficznym sensie”.

To znaczy konkretnie: bramka O3 dla `QW-2191` została domknięta według jawnych reguł i audytów w tym repo.

---

## 6. Mini-słownik licealisty

- **adversarial test** — test „złośliwy”: próbujemy celowo zepsuć wynik.  
- **residual ambiguity** — pozostała niejednoznaczność.  
- **theorem-level export** — formalny, jawny krok dowodowy, nie tylko wynik liczbowy.  
- **strict vs non-strict** — czy wynik wynika z rdzenia teorii, czy wymaga dodatkowego założenia.

---

## 7. Wniosek końcowy dla ucznia

W tej pracy zrobiono to, co robi dobra fizyka teoretyczna:

1. postawiono problem unikalności,
2. zbudowano mierzalne testy,
3. obalono słabe kandydaty,
4. utrzymano mocny kandydat,
5. domknięto formalne warunki,
6. sprawdzono niezależnie.

Dlatego w aktualnym stanie projektu można uczciwie powiedzieć:

\[
\boxed{\text{Bramka O3 dla } QW\text{-}2191 \text{ jest domknięta (strict, w rygorze pipeline).}}
\]

---

## 8. Co dalej?

Następny rozsądny krok po podręcznikowym domknięciu to nie „fajerwerki”, tylko higiena naukowa:

- regularne testy regresji,
- zewnętrzne replikacje,
- pilnowanie, żeby nowe zmiany nie odwróciły closure.

---

## 9. Historia wyprowadzenia domknięcia `QW-2191` do wersji 8 (prosty, ale ścisły skrót)

1. Najpierw rozdzielono porządek kerneli (legacy vs strict), żeby nie mieszać ról fizycznych bez jawnego mostu.
2. Potem zbudowano tor strict O3: katalog kandydatów, normalizacja, porównania pairwise, testy sweep/replay.
3. Po wykryciu słabych punktów (residual slot) przebudowano kandydaty (`v2`, `v3`, `v4`).
4. W wersji `v4` użyto jawnego score:

\[
s(\varphi,a)=(\varphi-0.40)+0.9(a-0.58)-0.6(\varphi-0.40)^2
\]

który decyduje o znaku gałęzi:

\[
\operatorname{sign}(s)=\begin{cases}+1,&s\ge0\\-1,&s<0\end{cases}
\]

5. Finalnie domknięto brak formalnego eksportu globalnego L2 i odświeżono checker (`5/5`).
6. Niezależna replikacja utrzymała wynik, więc w rygorze tej bramki status ustawiono na `CLOSED`.

## 10. Co kernel strict i score faktycznie wnoszą fizycznie?

Kernel strict
\[
K_{\text{strict}}(d)=\frac{\cos(\omega d+\phi)}{1+\beta d^{\eta}}
\]
porządkuje tłumienie i fazę oddziaływań w operacyjnej części modelu, a score `s(\varphi,a)` daje regułę selekcji znaku w mapie kandydatów.

Wkład fizyczny: redukcja niejednoznaczności selektora (unikalność gałęzi roboczej) przy zachowaniu ostrożności, że nie jest to jeszcze pełne „wszystko o SM+GR”.


---

## 11. Krótkie rozstrzygnięcie „tak/nie”

**TAK:** `QW-2191` jest domknięte w semantyce O3 Release 8.

**NIE:** to jeszcze nie jest pełne, „kernel-alone” rozwiązanie fundamentalnej niejednoznaczności z warstwy real Fourier basis i isotropy na pair planes.

Czyli uczciwy status jest podwójny: bramka O3 zamknięta, ale najmocniejsza wersja globalnego domknięcia nadal wymaga dodatkowych eksportów twierdzeń.


---

## 12. Domknięcie kernel-alone (P1340) — wersja uczniowska

Dodano pakiet `P1340`: domknięcie niejednoznaczności real Fourier basis + isotropy on pair planes zostało zrobione **warunkowo** przez jawne założenie selektora `KP1`.

Czyli: mamy formalne domknięcie kernel-alone w trybie „jeśli KP1, to unikalność”, a nie ukryte twierdzenie bez założeń.


---

## 13. Pełne domknięcie kernel-alone w ścieżce aksjomatycznej (P1342)

Podjęto decyzję: przyjmujemy jawny aksjomat `SB1` wyboru znaku i wtedy domykamy pełną niejednoznaczność kernel-alone.

To jest pełne domknięcie w ścieżce „z aksjomatem”, a nie ukryte twierdzenie „bez żadnych założeń”.


---

## 14. Pełne domknięcie strict przez źródło wewnętrzne (P1343)

Dodano `P1343`: model dostał wewnętrzne źródło selektora (`S_strict_internal_v1`), więc pełne domknięcie kernel-alone działa już także w wersji strict, nie tylko aksjomatycznej.


---

## 15. Walidacja obciążeniowa źródła strict (P1344)

Dodano `P1344`: źródło `S_strict_internal_v1` przeszło testy odporności, więc pełne domknięcie kernel-alone w wersji strict pozostaje utrzymane (z jawną polityką cofnięcia przy kontrprzykładzie).


---

## 16. Niezależna replikacja i próba obalenia (P1345)

Dodano `P1345`: wynik domknięcia został niezależnie powtórzony i nie znaleziono powtarzalnego kontrprzykładu, więc status strict pozostaje utrzymany.


---

## 17. Stabilność długohoryzontowa (P1346)

Dodano `P1346`: sprawdzono, że domknięcie strict utrzymuje się także w czasie i przy zmianach granic klasy dopuszczalnej (z nadal aktywną polityką cofnięcia).


---

## 18. Profesjonalny komunikat Release 8 (wersja publikacyjna)

Wersja 8 jest teraz opisana jako pełny pakiet naukowy dla klasy niejednoznaczności
(real Fourier basis + isotropy on pair planes) w ścieżce strict:

1. `P1343` — źródło selektora strict,
2. `P1344` — testy odporności + zasada cofnięcia,
3. `P1345` — niezależna replikacja i próba obalenia,
4. `P1346` — stabilność długohoryzontowa.

Status publikacyjny:

\[
\texttt{CLOSED\_FULL\_STRICT\_INTERNAL\_SOURCE\_V1\_LONG\_HORIZON\_STABLE}.
\]

To znaczy: dla tej klasy problemu teoria jest domknięta w rygorze strict opisanym w Release 8,
z zachowaniem jawnych granic i polityki cofnięcia przy przyszłym kontrdowodzie.


---

## 19. Dowiezienie globalnych bloków (P1347 + P1348)

Dodano dwa kluczowe pakiety:

1. `P1347` — formalna identyfikacja host-level w zadeklarowanym zakresie strict,
2. `P1348` — jedno globalne twierdzenie spinające cały łańcuch domknięcia.

W praktyce: globalne bloki z mapy blockerów zostały dowiezione dla zakresu Release 8.


---

## 20. Krok zewnętrzny: blind audit (P1349)

Dodano `P1349`: formalny protokół audytu zewnętrznego, który ma potwierdzić domknięcie przez niezależne zespoły i jawne kryteria pass/fail.


---

## 21. Aktualny stan (wersja jednoznaczna)

Na dziś dokumenty Release 8 należy czytać tak:

1. teoria jest domknięta w zadeklarowanym zakresie strict Release 8,
2. istnieje formalne globalne twierdzenie domknięcia (`P1348`),
3. jedyny krok „po domknięciu” to wykonanie zewnętrznego blind audytu (`P1349`).

Czyli: domknięcie jest pełne w obecnym zakresie wiedzy, a teraz czekamy na zewnętrzne potwierdzenie.
