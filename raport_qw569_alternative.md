# Raport QW-569: Alternatywne Metody Pomiaru

**Data:** 2025-12-05  
**Status:** CZĘŚCIOWY SUKCES / PORAŻKA  
**Cel:** Sprawdzenie czy alternatywne metody (efekty) działają lepiej niż bezpośrednie pola na N=400.

---

## Wyniki Testów

### **1. Propagation Delay (v(r))**
- **Metoda:** Pomiar czasu przylotu impulsu falowego.
- **Wynik:** **Niewystarczająca ilość danych** (1 punkt).
- **Przyczyna:** Impuls zanikł lub rozproszył się zbyt szybko w szumie sieci N=400.

### **2. Correlation Functions (Time Dilation)**
- **Metoda:** Pomiar czasu zaniku korelacji $\langle \psi_A(0) | \psi_A(t) \rangle$.
- **Wynik:** $\tau_A = 0.100$, $\tau_B = 0.100$.
- **Ratio:** $\gamma = 1.0000$ (Teoria: 1.0043).
- **Błąd:** 100% (brak wykrycia efektu).
- **Wniosek:** Dekoherecja jest natychmiastowa ($\tau=0.1$ to 1 krok symulacji!). Szum sieci dominuje nad subtelnym efektem dylatacji.

### **3. Angular Momentum (Frame Dragging)**
- **Metoda:** Pomiar $\langle L_z \rangle$ w powłokach.
- **Wynik:** $n=0.000$ (Teoria: $n=2.0$).
- **Wniosek:** Brak sygnału wleczenia układu. Wir zanika zbyt szybko lub jest zbyt słaby na N=400.

---

## Diagnoza

Mimo zmiany metody na "pomiar efektów", sieć N=400 okazała się **zbyt mała i zbyt szumna** dla tych subtelnych zjawisk:
1.  **Szybka dekoherencja:** Stany kwantowe tracą koherencję w 1-2 kroki, uniemożliwiając pomiar dylatacji czy propagacji.
2.  **Słaby sygnał:** Efekty rzędu 0.4% giną w szumie dyskretyzacji.

**Potwierdzenie:**
Hipoteza "Effects > Fields" jest prawdziwa dla **dużych efektów** (orbity, QW-565), ale dla **mikro-efektów** (dylatacja, dragging) potrzebna jest **większa rozdzielczość** (mniejszy szum).

---

## Rekomendacja: QW-570

Zgodnie z sugestią użytkownika, należy przejść do **QW-570: Frame Dragging Enhanced**:
1.  **N=1500+:** Zmniejszenie szumu dyskretyzacji.
2.  **Dłuższa ewolucja (500+ kroków):** Pozwolenie na rozwinięcie się efektów.
3.  **Lepsza inicjalizacja wiru:** Aby przetrwał dekoherencję.
4.  **Metoda Momentu Pędu:** Użycie metody z QW-569 (jest lepsza niż gradient fazy), ale na większej sieci.

---

**Podsumowanie:**
QW-569 pokazał, że sama zmiana metody nie wystarczy na N=400. Skala sieci MA znaczenie dla mikro-efektów (zmniejszenie szumu). Przechodzimy do QW-570.
