# RAPORT NAUKOWY: REBOOT BADAŃ (QW-735 - QW-754)

> [!IMPORTANT]
> **WERDYKT: CHAOS KWANTOWY I WYMIAR 2.3D**
> Pełna symulacja numeryczna ($N=2000$) obala hipotezę o "1D Włóknach". W reżimie połączonym (Percolation) Eter zachowuje się jak **Chaotyczna Piana** o wymiarze spektralnym $d_s \approx 2.27$.

## 1. Rygor Metodologiczny ("Zero Assumption")
Dane pochodzą wyłącznie z bezpośredniej diagonalizacji Laplasjanu i śledzenia geodezyjnych na grafie o 1989 węzłach (Giant Component).
*   **Engine:** `scipy.sparse.linalg.eigsh` (Exact Diagonalization).
*   **Metric:** Dijkstra Geodesics.

## 2. Kluczowe Wyniki (The Hard Data)

### A. Wymiar Rzeczywistości (QW-736, QW-741)
*   **Spectral Dimension:** $d_s \approx 2.27$.
*   **Box Dimension:** $D_0 \approx 2.08$.
*   **Interpretacja:** Przestrzeń jest efektywnie **2-wymiarowa** (jak pognieciona kartka papieru), a nie 3-wymiarowa ani 1-wymiarowa. To sugeruje model **Grawitacji Kwantowej Liouville'a** (2D Quantum Gravity).

### B. Chaos Kwantowy (QW-744)
*   **Level Spacing Ratio ($r$):** $0.511$.
*   **Teoria:**
    *   Poisson (Integrowalny/Regularny): $r \approx 0.38$.
    *   GOE (Chaotyczny/Wigner-Dyson): $r \approx 0.53$.
*   **Werdykt:** Eter jest systemem **Kwantowo Chaotycznym** (blisko GOE). Informacja w nim ulega natychmiastowemu "wymieszaniu" (scrambling). To jest cecha Czarnych Dziur!

### C. Anizotropia i Prędkość (QW-739, QW-743)
*   **Anisotropy Index:** $0.136$ (Niska).
*   **Signal Speed:** $c \approx 0.62$ (geometryczna).
*   **Wniosek:** Mimo fraktalnej struktury, w skali "dużej" ($N=2000$) przestrzeń staje się izotropowa. Nie ma "uprzywilejowanego kierunku".

## 3. Co Obalono? (Failures)
1.  **Hipoteza 1D Włókien (Poprzednia):** Obalona. Przy $R=1.0$ (spójny wszechświat) wymiar skacze do >2. Włókna były artefaktem zbyt małego promienia łączenia ($R=0.8$).
2.  **Pusta Przestrzeń 3D:** Obalona. Wymiar $2.27 \ll 3$. Eter jest bardzo "porowaty" lub "płaski".

## 4. Konkluzja
Przeszliśmy z modelu "Sznurka" do modelu **"Kwantowej Membrany" (2D Gravity)**, która wykazuje cechy chaosu kwantowego.
Następny krok: Sprawdzić, czy na tej membranie mogą istnieć stabilne wzbudzenia (Anyony?).
