# RAPORT NAUKOWY: DYNAMIKA WĘZŁÓW I NIESZCZĘSNA PRAWDA (QW-787 - QW-806)

> [!WARNING]
> **WERDYKT: NIESTABILNOŚĆ**. Hipoteza, że "cząstki to trwałe węzły w sieci", napotkała zderzenie z rzeczywistością termiczną.

## 1. Metodologia (Rygor Naukowy)
1.  **Detekcja Cykli:** Wykorzystano algorytm oparty na Drzewie Rozpinającym (MST) do znalezienia fundamentalnych pętli bazy homologii ($H_1$).
2.  **Liczba Splątania Gaussa ($Lk$):** Obliczona numerycznie jako podwójna całka dyskretna po krawędziach grafu.
    $$ Lk = \frac{1}{4\pi} \sum_{i,j} \frac{\mathbf{r}_{ij} \cdot (\Delta \mathbf{r}_i \times \Delta \mathbf{r}_j)}{|\mathbf{r}_{ij}|^3} $$
3.  **Symulacja Monte Carlo:** Wprowadzono szum termiczny (przełączanie krawędzi) i mierzono czas zaniku topologii.

## 2. Wyniki "The Ugly Truth" (Brzydka Prawda)

### A. Cząstki nie są Trwałe (QW-796)
*   **Lifetime:** ~120 kroków symulacji.
*   **Wniosek:** Węzeł w Eterze Nadsolitona nie jest jak "proton" (trwały wiecznie). Jest jak "wir w wodzie" – powstaje i znika.
*   **Problem:** Jak zbudować stabilną materię z nietrwałych fluktuacji?

### B. Kwantyzacja jest Przybliżona (QW-801)
*   **Błąd Kwantyzacji:** $\approx 0.096$ (ok. 10%).
*   **Relacja:** $Lk \approx 0.18$ (max).
*   **Wniosek:** Węzły są "słabe". Rzadko zdarza się pełne, twarde splątanie ($Lk=1$). To raczej "zahaczenia" niż węzły żeglarskie.

### C. Brak Pauli'ego (QW-798)
*   **Status:** Nie wykryto topologicznego zakazu nakładania się pętli (Exclusion = None).
*   **Problem:** To sugeruje, że Eter jest naturalnie bozonowy. Fermiony (elektrony) muszą być czymś znacznie bardziej złożonym (np. twistem, a nie pętlą).

## 3. Co to oznacza dla Teorii?
Musimy zaakceptować **Paradygmat Fluktuacyjny**.
Materia, którą widzimy, nie jest sumą statycznych węzłów. Jest **dynamicznym stanem równowagi** procesów wiązania i rozwiązywania.
*   Cząstka to "Soliton Topologiczny" – fala, która *odnawia* swoją strukturę węzła szybciej, niż szum go niszczy.
*   Stabilność protonu wymaga **mechanizmu aktywnego podtrzymywania** (Rezonans), a nie tylko pasywnej topologii.

**Werdykt:** Topologia sama w sobie nie wystarcza do stabilizacji materii. Potrzebna jest "Termodynamika Nierównowagowa" (Non-equilibrium Thermodynamics) lub Rezonans (QW-747+).
