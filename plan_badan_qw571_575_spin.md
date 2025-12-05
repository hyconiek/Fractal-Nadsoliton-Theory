# Plan Badań QW-571 do QW-575: Spin Networks & Quantum Gravity

**Cel:** Przejście od modelu skalarnego ($\Psi$) do modelu wektorowego/tensorowego (Spin Networks), aby naprawić brak rotacji (Frame Dragging) i połączyć teorię FIN z Pętlową Grawitacją Kwantową (LQG).

**Uzasadnienie:**
Porażka QW-570 (Frame Dragging) wykazała, że skalarna sieć neuronowa nie posiada stopni swobody niezbędnych do przenoszenia momentu pędu. Potrzebujemy "kwantowej geometrii" opartej na spinach.

---

## Lista Eksperymentów

### **QW-571: Spin Network Initialization**
*   **Cel:** Zdefiniowanie sieci, w której każdy węzeł posiada wektor spinu $\vec{S}$ (lub spinor SU(2)), a krawędzie niosą reprezentacje nieprzywiedlne (kolorowanie krawędzi).
*   **Metoda:** Rozszerzenie klasy `Node` o obiekt `SpinState`.
*   **Oczekiwany Wynik:** Stabilna inicjalizacja sieci spinowej.

### **QW-572: Spin-Spin Interaction (Heisenberg Model)**
*   **Cel:** Wprowadzenie dynamiki opartej na sprzężeniu spinów: $H = -J \sum \vec{S}_i \cdot \vec{S}_j$.
*   **Pytanie:** Czy w sieci fraktalnej wyłania się uporządkowanie magnetyczne (ferro/antyferro) czy ciecz spinowa (Quantum Spin Liquid)?
*   **Znaczenie:** Ciecz spinowa jest kandydatem na "próżnię kwantową".

### **QW-573: Emergent Geometry (Area Operator)**
*   **Cel:** Sprawdzenie, czy geometria (powierzchnia/objętość) wyłania się ze spinów, jak w LQG.
*   **Metoda:** Obliczenie widma operatora pola powierzchni: $A \propto \sum \sqrt{j(j+1)}$.
*   **Oczekiwany Wynik:** Dyskretne widmo powierzchni (kwantyzacja przestrzeni).

### **QW-574: Frame Dragging with Spin (Re-test QW-570)**
*   **Cel:** Powtórzenie testu QW-570 w modelu spinowym.
*   **Metoda:** Rotująca masa wymusza precesję spinów w otoczeniu.
*   **Oczekiwany Wynik:** **SUKCES**. Przeniesienie momentu pędu przez sprzężenie spin-spin (spin current).

### **QW-575: Quantum Graphity (Dynamic Rewiring)**
*   **Cel:** Symulacja ewolucji samej topologii sieci (dodawanie/usuwanie krawędzi) pod wpływem stanu spinowego.
*   **Hipoteza:** "Geometria to stan niskoenergetyczny sieci spinowej". W wysokich temperaturach przestrzeń znika (faza pre-geometryczna).

---

## Metodologia
*   **Język:** Python + NumPy (macierze Pauliego).
*   **Skala:** N=500-1000 (obliczenia na spinach są kosztowne).
*   **Weryfikacja:** Porównanie z przewidywaniami LQG (widmo powierzchni) i OTW (frame dragging).
