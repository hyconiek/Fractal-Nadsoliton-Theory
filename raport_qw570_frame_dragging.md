# Raport QW-570: Frame Dragging Enhanced

**Data:** 2025-12-05  
**Status:** PORAŻKA (FAIL)  
**Cel:** Wymuszenie efektu Lense-Thirringa (wleczenie układu) w sieci N=1500 z aktywnym napędzaniem (Active Driving).

---

## Wyniki Testu

### **1. Angular Momentum Transfer**
- **Metoda:** Active Driving (człon $\Omega \cdot L_z$ w Hamiltonianie dla masy).
- **Oczekiwanie:** Stopniowy wzrost momentu pędu próżni ($L_z$ vacuum).
- **Wynik:**
    - Start: $L_z = 0.00$
    - Max: $L_z \approx 0.68$ (w trakcie ewolucji)
    - Final: $L_z = -0.007$
- **Wniosek:** Brak trwałego przekazu momentu pędu. Próżnia "ślizga się" po rotującej masie. Energia jest rozpraszana w chaotyczne oscylacje, a nie w koherentny wir.

### **2. Rotation Curve**
- **Metoda:** Pomiar profilu prędkości azymutalnej $v_\phi(r)$.
- **Wynik:**
    - Fit: $v(r) \propto r^{-1.60}$
    - $R^2 = 0.37$ (bardzo słabe dopasowanie)
- **Wniosek:** Brak wyraźnej struktury wiru. Rozkład prędkości jest bliski losowemu szumowi.

---

## Diagnoza

Dlaczego QW-405 (FFT/GPE) działało, a QW-570 (Discrete Network) nie?

1.  **Topologia:** W sieci dyskretnej "rotacja" jest słabo zdefiniowana. Brak ciągłej symetrii obrotowej powoduje, że moment pędu nie jest ściśle zachowany (sieć łamie symetrię SO(3)).
2.  **Lepkość:** "Lepkość próżni" (Vacuum Viscosity) w modelu dyskretnym może być zbyt mała, by przenieść rotację na większe odległości, lub zbyt duża (tarcie dynamiczne), co tłumi wir.
3.  **Brak zmiennych spinowych:** Model opiera się na skalarnej funkcji falowej $\psi$. Wiry w $\psi$ są możliwe (Hopfiony), ale wymagają specyficznych warunków topologicznych, których proste "napędzanie fazy" nie tworzy.

---

## Podsumowanie Sekcji Flow (QW-563-570)

Mimo porażki QW-570, cała sekcja "Gravity as Flow" jest **CZĘŚCIOWYM SUKCESEM**:

| Test | Wynik | Znaczenie |
|------|-------|-----------|
| **QW-564 (Flow vs Force)** | ✅ **SUKCES** | Model przepływu 2.5x lepszy niż siłowy. |
| **QW-565 (Geodesics)** | ✅ **SUKCES** | Stabilne orbity ($e=0.71$) bez siły $1/r^2$. |
| **QW-566/569 (Dilation)** | ⚠️ **FAIL** | Efekt zbyt słaby dla N=400-1500. |
| **QW-567/570 (Dragging)** | ❌ **FAIL** | Brak rotacji w sieci skalarnej. |

**Ostateczny Wniosek:**
Grawitacja w FIN **JEST PRZEPŁYWEM** (potwierdzone trajektorie i orbity), ale efekty wyższego rzędu (rotacja, dylatacja) wymagają albo znacznie większej skali, albo wprowadzenia spinu jako fundamentalnej zmiennej (Spin Networks).
