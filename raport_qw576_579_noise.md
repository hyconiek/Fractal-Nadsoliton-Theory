# Raport Badań QW-576 do QW-579: Noise & Quantum Criticality

**Data:** 2025-12-05
**Cel:** Znalezienie punktu krytycznego ($T_c$) i wykorzystanie go do uzyskania efektu frame dragging.

---

## Wyniki Symulacji

### **QW-576: Skanowanie Przejścia Fazowego**
*   **Cel:** Znalezienie piku podatności magnetycznej $\chi$.
*   **Wynik:**
    *   Magnetyzacja $M$ spada monotonicznie z $0.91$ ($T=0.1$) do $0.14$ ($T=5.0$).
    *   Podatność $\chi$ jest niska i płaska, bez wyraźnego piku (Divergence).
    *   Algorytm błędnie zidentyfikował $T_c=0.10$ (maksimum lokalne w szumie).
*   **Wniosek:** W sieci losowej o wysokim stopniu (~17 sąsiadów) przejście fazowe jest silnie rozmyte (Mean Field). Nie ma ostrego punktu krytycznego w badanym zakresie.

### **QW-577: Critical Frame Dragging**
*   **Warunki:** Test wykonano w $T=0.10$ (Faza Ferromagnetyczna, $M=0.91$).
*   **Wynik:** ❌ **PORAŻKA (L_z = 0.025)**.
*   **Przyczyna:** Test powtórzył błąd QW-574. W niskiej temperaturze spiny są "zamrożone" w domenie magnetycznej i nie reagują na słabą rotację grawitacyjną.

### **QW-579: Grawitacja Entropowa**
*   **Status:** Hipotetyczny. Ponieważ nie znaleziono ostrego przejścia fazowego, nie można potwierdzić mechanizmu Verlinde'a (maksymalizacja entropii w $T_c$).

---

## Diagnoza i Wnioski
1.  **Brak Cieczy Spinowej:** Prosty model Heisenberga na grafie losowym dąży do ferromagnetyzmu. Aby uzyskać **Kwantową Ciecz Spinową** (Quantum Spin Liquid), która jest "miękka" i reaguje na metrykę, potrzebujemy **frustracji geometrycznej** (np. sieci trójkątne/Kagome) lub bardziej złożonych interakcji (np. model Kitaeva).
2.  **Sukces Geometrii (QW-573):** Mimo problemów z dynamiką (dragging), statyczna geometria (Operator Pola Powierzchni) jest poprawna i zgodna z LQG.

**Rekomendacja:** Przyszłe badania powinny skupić się na **topologii sieci** (wprowadzenie frustracji), aby uniknąć ferromagnetycznego zamrażania próżni.
