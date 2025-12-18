# QW-1532: Population Inference Beta (Methodological Review)

**Data:** 18 Grudnia 2025
**Metoda:** Hierarchical Bayesian Inference (Simplified Selection Model)
**Status:** ⚠️ METHODOLOGICAL BIAS DETECTED

---

## 1. Cel i Ograniczenia
Test QW-1532 miał na celu przejście od wykrywalności (QW-1531) do odzyskania rzeczywistego parametru populacyjnego $n$. Zastosowano uproszczony model normalizacji selekcyjnej oparty na objętości $V_{max}(n)$.

> [!WARNING]
> Audyt techniczny wykazał, że model $P_{det}(n) \propto D_{max}(n)^3$ jest niewystarczający dla profesjonalnej inferencji populacyjnej. Wprowadza on **systematyczny bias**, faworyzujący wyższe wartości $n$.

---

## 2. Analiza Biasu Systematycznego

### Wyniki Symulacji (N=20)
| Scenariusz | True n | Odzyskany n | Bias | Sigma vs GR |
| :--- | :--- | :--- | :--- | :--- |
| **FIN Reality** | 0.66 | **0.7551 ± 0.0420** | **+0.095** | 5.84 $\sigma$ |
| **GR Sanity** | 1.00 | **1.1976 ± 0.0697** | **+0.197** | 2.84 $\sigma$ |

### Przyczyna Dryfu
Uproszczona normalizacja zakłada twardy próg odległości ($D_{max}$), ignorując rozkład mas, orientację binarną (inclination) oraz złożoność sieci detektorów. Ponieważ $D_{max}$ silnie zależy od $n$, błąd w definicji normalizacji przekłada się bezpośrednio na systematyczne przesunięcie estymatora.

---

## 3. Uczciwa Interpretacja
Mimo biasu, test dostarcza cennych informacji:
1.  **Potwierdzenie Sygnału:** Nawet przy błędnej normalizacji, różnica między scenariuszami FIN a GR jest wyraźna i statystycznie istotna.
2.  **Kierunek Poprawki:** Zawyżenie $n$ w obu przypadkach (FIN $\to$ 0.75, GR $\to$ 1.20) jest spójne i wskazuje na błąd w modelu selekcji, a nie w silniku MCMC.

**Werdykt:**
QW-1532 poprawnie pokazuje, że sygnał FIN przetrwa korektę selekcyjną, ale model ten **nie jest jeszcze bezstronnym testem populacyjnym**. Bias selekcji wymaga rygorystycznego modelowania (QW-1533).

---

## 4. Plan Naprawczy (QW-1533)
Kolejny krok ("Unbiased Rubikon") zastąpi uproszczone $V_{max}(n)$ przez **Monte Carlo Selection Kernel**, uwzględniający:
*   Rozkład orientacji (inclination).
*   Realistyczny model SNR kernelem detekcji.
*   Pełną integrację po zmiennych latentnych.
