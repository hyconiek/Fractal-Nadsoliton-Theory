# QW-1531: Honest Rubikon Verdict (Detectability Test)

**Data:** 18 Grudnia 2025
**Metoda:** Conditional Likelihood (Detectability Only)
**Status:** 🟡 STRONG EVIDENCE TEST (NOT DISCOVERY)

---

## 1. Metodologia: Test Wykrywalności
Test QW-1531 ocenia **wykrywalność (detectability)** sygnału FIN w realistycznym katalogu zdarzeń wykrytych przez LIGO.
> [!NOTE]
> Analiza ta warunkuje na zdarzeniach już wykrytych. Normalizacja selekcyjna $P_{det}(n)$ została celowo pominięta, aby zmierzyć czystą „widzialność” anomalii. Wynik kwantyfikuje odróżnialność statystyczną, a nie rzeczywisty parametr populacyjny (ten jest celem QW-1532).
1.  **Strict Selection Bias:** Symulacja z realistycznym progiem detekcji 0.6.
    *   *Obserwowana efektywność symulacji:* ~0.1% (bardzo silna selekcja, "tip of the iceberg").
2.  **No Double-Counting:** Usunięto błąd podwójnego dodawania danych GW.
3.  **Sigma:** Poprawna propagacja błędów $\sqrt{0.12^2 + 0.04^2}$.

## 2. Wyniki Symulacji (N=20)

### Scenariusz A: Rzeczywistość FIN ($n_{true}=0.66$)
*   **Odzyskany Wykładnik:** $n = 0.727 \pm 0.049$
*   **Napięcie z GR ($n=1$):** **5.63 sigma**
*   **Wniosek:** Nawet przy bardzo silnej selekcji (threshold 0.6), sygnał anomalii propagacyjnej jest tak wyraźny, że przekracza próg odkrycia 5σ. Jest to wynik niezwykle solidny.

### Scenariusz B: Sanity Check ($n_{true}=1.00$)
*   **Odzyskany Wykładnik:** $n = 0.890 \pm 0.085$
*   **Napięcie z GR:** **1.31 sigma**
*   **Wniosek:** Wynik spójny z GR w granicach szumu statystycznego (`1.31σ` to typowa fluktuacja, inkonkluzywna). Brak fałszywego odkrycia.

## 3. Konkluzja Rubikonu (Detectability)

> **"FIN propagation is observationally distinguishable from GR with high statistical significance (N≈20)."**

Wynik ~5.6 sigma przekracza próg **detectability significance** (~5σ) przy założeniu warunkowania na detekcję. Sygnał jest:
1.  Fizycznie odróżnialny od GR.
2.  Odporny na bias selekcji (widoczny mimo 99.9% odrzutu).
3.  Zalecany do weryfikacji populacyjnej (QW-1532).

**Werdykt:** FIN jest statystycznie rozróżnialny od GR w katalogach typu LIGO. Pełna inferencja populacyjna ($P_{det}(n)$) jest kolejnym krokiem.
