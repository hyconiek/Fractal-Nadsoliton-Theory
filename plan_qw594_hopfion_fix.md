# Plan Badań QW-594: Revised Hopfion Stability (Attractor Dynamics)
**Data:** 2025-12-05
**Cel:** Ponowna próba weryfikacji H4 (Cząstki = Wiry) po porażce QW-590.

---

## 1. Analiza Porażki QW-590
*   **Objawy:** Wartości `NaN` w symulacji, natychmiastowy rozpad topologii.
*   **Przyczyny:**
    1.  **Osobliwości:** Obliczanie ładunku topologicznego dzieliło przez `rho` (gęstość), która w rdzeniu wiru dąży do zera.
    2.  **Niestabilność Numeryczna:** Brak odpowiedniej regularyzacji w ewolucji (zbyt duże DT lub brak dyfuzji stabilizującej).
    3.  **Złe Założenia Dynamiki:** Użyto tylko potencjału efektywnego. Nowa wiedza (QW-558) mówi, że Nadsoliton to **Atraktor** ($dA/dt = A - A^3$).

## 2. Nowe Podejście (QW-594)

### A. Dynamika Atraktora
Zamiast równania Schrödingera z potencjałem, użyjemy **Ginzburg-Landau-like Attractor Equation**:
$$\frac{\partial \psi}{\partial t} = (1 - |\psi|^2)\psi + i \nabla^2 \psi + i V_{eff} \psi$$
*   Człon $(1 - |\psi|^2)\psi$: Wymusza relaksację amplitudy do 1 (próżnia), ale pozwala na zero w rdzeniu wiru. To jest mechanizm stabilizujący (Negative Feedback).
*   Człon $i \nabla^2 \psi$: Kinetyka / Dyspersja.
*   $V_{eff}$: Wpływ warstw fraktalnych (Back-reaction).

### B. Numeryka
1.  **Regularyzacja Ładunku Topologicznego:**
    Zamiast dzielić przez `rho`, użyć formy całkowej, która jest bezpieczna numerycznie, lub dodać `epsilon` do mianownika.
2.  **Mniejszy krok czasowy:** `DT = 0.01` (zamiast 0.05).
3.  **Soft-Core Initialization:** Inicjalizacja z łagodniejszym profilem, aby uniknąć "szoku" na początku symulacji.

### C. Hipoteza Badawcza
Jeśli H4 jest prawdziwa, to **dynamika atraktora (Network Self-Optimization)** powinna "naprawić" drobne błędy i utrzymać topologię hoffiona, podczas gdy czysta mechanika kwantowa (unitary) pozwala na chaotyczny rozpad.

## 3. Plan Implementacji
1.  Stworzyć `QW-594_Hopfion_Revised.py`.
2.  Zaimplementować równanie Ginzburga-Landaua.
3.  Zmierzyć nie tylko Q (ładunek), ale i **Energię** układu (czy maleje do minimum lokalnego?).

---
