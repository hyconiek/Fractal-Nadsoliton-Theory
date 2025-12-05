# Raport Badań QW-515 do QW-519: Effective Fractal Coupling

**Data:** 2025-12-04
**Cel:** Szybka weryfikacja numeryczna (bez fittingu) hipotez o "Efektywnym Sprzężeniu Fraktalnym" przy zamrożonych parametrach.

---

## Wyniki Badań

### **QW-515: Test Odwrotnej Hierarchii (Echo)**
*   **Cel:** Sprawdzić, czy `K(d)` rośnie dla dużych `d` (odbicie siły).
*   **Wynik:** ✅ **SUKCES**.
    *   `|K(6)| = 1.31`
    *   `|K(12)| = 2.14`
    *   **Wniosek:** Istnieje silne "echo" na 12. oktawie (`|K(12)| > |K(6)|`). Potwierdza to mechanizm bezpośredniego sprzężenia góry z dołem (Inverse Hierarchy).

### **QW-516: Efektywny Potencjał Wodoru (N=10)**
*   **Cel:** Sprawdzić, czy potencjał `-K(r) - K(10)/r` generuje widmo $1/n^2$.
*   **Wynik:** 🔴 **PORAŻKA**.
    *   Stosunek $E_2/E_1 = 0.96$ (Oczekiwane: $0.25$).
    *   **Wniosek:** Prosty model potencjału efektywnego z tłem od 10. warstwy nie odtwarza struktury atomu wodoru. Potrzebna jest pełna symulacja rezonansu sieci (jak w planie QW-500).

### **QW-517: Masa z Tłumienia Skalowego**
*   **Cel:** Sprawdzić, czy tłumienie przez 10 warstw ($\prod (1+\beta i)^{-1}$) daje stosunek masy protonu do Plancka ($10^{-19}$).
*   **Wynik:** 🔴 **PORAŻKA (Ilościowa)**.
    *   Tłumienie dla $N=10$ wynosi $\sim 0.59$.
    *   Wymagane tłumienie $10^{-19}$ osiągane jest dopiero dla **$N=108$**.
    *   **Wniosek:** Hipoteza "12 oktaw" jest niewystarczająca do wyjaśnienia masy protonu przez proste tłumienie. Struktura musi być głębsza ($N \approx 100$) lub tłumienie silniejsze (nieliniowe).

### **QW-518: Stała Hubble'a z Efektu Casimira**
*   **Cel:** Powiązać $H$ z gęstością energii próżni dla $N=30$.
*   **Wynik:** 🔴 **PORAŻKA**.
    *   Obliczone $H \sim 10^{-90}$ (jednostki Plancka).
    *   Obserwowane $H_0 \sim 10^{-61}$.
    *   Rozbieżność rzędu 30 rzędów wielkości.
    *   **Wniosek:** Model liniowego sumowania energii próżni ($\sum K$) jest błędny.

### **QW-519: Niezmienniczość Stałej Struktury ($\alpha$)**
*   **Cel:** Sprawdzić zachowanie $\alpha$ przy skalowaniu $\beta$.
*   **Wynik:** 🟡 **CZĘŚCIOWE / INSIGHT**.
    *   Jeśli skalujemy tylko $\beta \to \beta/S$, to $\alpha^{-1}$ rośnie do nieskończoności.
    *   Jeśli skalujemy również $\alpha_{geo} \to \alpha_{geo}/S$, to $\alpha^{-1}$ pozostaje stabilne ($\approx 138.6$).
    *   **Wniosek:** Aby stała struktury subtelnej była niezmiennikiem fraktalnym, parametr geometryczny $\alpha_{geo}$ musi być skalowalny (zależeć od skali), a nie stały.

---

## Podsumowanie
Badania wykazały, że **proste modele analityczne (QW-516, 517, 518) zawodzą**. Jedynie **QW-515 (Echo)** potwierdziło przewidywania teorii.
Kluczowym odkryciem z **QW-519** jest konieczność skalowania parametru $\alpha_{geo}$ wraz z $\beta_{tors}$, aby zachować spójność fizyczną (stałą $\alpha_{EM}$).

**Rekomendacja:** Porzucić proste modele "jednowarstwowe" na rzecz pełnych symulacji sieciowych (zgodnie z planem Phase XVII).
