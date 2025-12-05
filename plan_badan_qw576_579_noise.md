# Plan Badań QW-576 do QW-579: Noise & Quantum Criticality

**Cel:** "Rozmrożenie" próżni spinowej (która w QW-572 okazała się sztywnym ferromagnetykiem) poprzez wprowadzenie kontrolowanego szumu (temperatury), aby osiągnąć stan **Kwantowej Cieczy Spinowej** (Quantum Spin Liquid) lub **Punktu Krytycznego**.

**Uzasadnienie:**
W QW-574 frame dragging nie zadziałał ($L_z \approx 0.05$), bo spiny były zablokowane ($M=1.0$). W punkcie krytycznym podatność magnetyczna $\chi$ dąży do nieskończoności, co powinno umożliwić nawet słabemu polu grawitacyjnemu (rotacji) pociągnięcie spinów.

---

## Lista Eksperymentów

### **QW-576: Phase Transition Scan (Finding T_c)**
*   **Cel:** Znalezienie temperatury krytycznej $T_c$, w której następuje przejście z fazy uporządkowanej (Ferro) do nieuporządkowanej (Para).
*   **Metoda:** Symulacja Monte Carlo dla różnych $T$. Pomiar magnetyzacji $M$ i podatności $\chi = (\langle M^2 \rangle - \langle M \rangle^2) / T$.
*   **Oczekiwany Wynik:** Pik podatności w $T_c$. To jest nasz "Sweet Spot".

### **QW-577: Critical Frame Dragging (The Fix)**
*   **Cel:** Powtórzenie testu QW-574 (Frame Dragging) dokładnie w temperaturze $T \approx T_c$.
*   **Hipoteza:** W punkcie krytycznym korelacje są dalekozasięgowe, a sztywność znika. Rotacja powinna wywołać silną reakcję spinów.
*   **Oczekiwany Wynik:** Znaczący wzrost $L_z$ (np. $> 0.2$).

### **QW-578: Stochastic Resonance (Noise Benefit)**
*   **Cel:** Sprawdzenie, czy szum pomaga w stabilizacji struktur (Nadsolitona).
*   **Metoda:** Dodanie szumu do ewolucji Nadsolitona (QW-558) i pomiar czasu życia/stabilności.
*   **Hipoteza:** Umiarkowany szum (Stochastic Resonance) zapobiega utknięciu w lokalnych minimach i wzmacnia sygnał.

### **QW-579: Emergent Gravity from Noise (Entropic Force)**
*   **Cel:** Czy grawitacja jest siłą entropową wynikającą z szumu?
*   **Metoda:** Pomiar siły przyciągania w funkcji temperatury.
*   **Hipoteza:** Grawitacja Verlinde'a wymaga $\Delta S > 0$. W $T=0$ (Ferro) entropia jest zerowa. W $T_c$ entropia jest maksymalna. Grawitacja powinna być najsilniejsza w pobliżu przejścia fazowego.

---

## Metodologia
*   **Model:** Spin Network (z QW-571) + Termostat (Langevin Dynamics lub Metropolis).
*   **Parametr Sterujący:** Temperatura $T$ (reprezentująca szum kwantowy/fluktuacje próżni).
