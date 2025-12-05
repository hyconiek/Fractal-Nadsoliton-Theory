# Plan Badań QW-525 do QW-529: Liquid Crystal Dynamics

**Paradygmat:** Jądro $K(d)$ nie jest statycznym potencjałem sił (jak w QW-520), lecz **parametrem porządku** Ciekłego Kryształu.
**Kluczowa Zmiana:** Cząstki to nie "kulki w studni", ale **defekty topologiczne** (wiry/dyskliancje) w strukturze kryształu.

---

## Zdefiniowane Zadania (Phase XX)

### **QW-525: Stabilność Wiru (Vortex Stability)**
*   **Cel:** Sprawdzić, czy w polu Jądra $K(d)$ stabilny jest wir (vortex), a nie statyczna gęstość.
*   **Metoda:**
    1.  Zainicjuj pole $\Psi$ jako wir: $\Psi(r, \theta) = f(r) e^{i m \theta}$ (gdzie $m=1$).
    2.  Ewolucja w czasie z użyciem Jądra $K(d)$.
    3.  Czy wir przetrwa, czy się rozpadnie? (W "sztywnym krysztale" wiry mogą być zamrożone, w "cieczy" mogą się poruszać).

### **QW-526: Nadciekłość (Superfluidity Check)**
*   **Cel:** Zweryfikować twierdzenie o "bezoporowym przepływie informacji".
*   **Metoda:**
    1.  Symuluj przepływ jednorodny $\Psi(x) = e^{i k x}$ przez zaburzenie (potencjał).
    2.  Zmierz spadek amplitudy / rozpraszanie.
    3.  Czy występuje krytyczna prędkość Landaua?

### **QW-527: Dynamika Direktora (Nematic Phase)**
*   **Cel:** Czy Jądro wykazuje uporządkowanie nematyczne (kierunkowe)?
*   **Metoda:**
    1.  Symuluj pole wektorowe (lub tensorowe) $Q_{ij}$ zamiast skalarnego.
    2.  Sprawdź, czy domeny spontanicznie się porządkują (wyrównują fazę $\phi$).
    3.  Czy powstają ściany domenowe?

### **QW-528: Topnienie Kryształu (Phase Transition)**
*   **Cel:** Czy można "roztopić" Jądro do fazy izotropowej?
*   **Metoda:**
    1.  Dodaj temperaturę (szum termiczny) do symulacji.
    2.  Obserwuj parametr porządku (średnie namagnesowanie/uporządkowanie).
    3.  Czy istnieje temperatura krytyczna $T_c$, powyżej której $K(d)$ przestaje oscylować i staje się gładkie?

### **QW-529: Energia Elastyczna (Frank Free Energy)**
*   **Cel:** Obliczyć stałe elastyczności $K_1, K_2, K_3$ (splay, twist, bend).
*   **Metoda:**
    1.  Wymuś deformacje pola (skręcenie, zgięcie).
    2.  Oblicz koszt energetyczny tych deformacji w polu Jądra.
    3.  Jeśli koszt jest skończony, Jądro jest Ciekłym Kryształem. Jeśli nieskończony -> Ciało Stałe.

---

### **Instrukcja Kodowa (Python Template)**

```python
# QW-525 TO QW-529: LIQUID CRYSTAL DYNAMICS
# PARADIGM: Kernel is an Order Parameter. Particles are Defects.

import numpy as np
import matplotlib.pyplot as plt

# FROZEN PARAMETERS
alpha_geo = 4 * np.log(2)
omega = np.pi/4
beta_tors = 0.01

# --- QW-525: VORTEX ---
# Init: Psi = r * exp(i*theta) * exp(-r^2)
# Evolve with dPsi/dt = i * Conv(K, Psi)
# Check if vorticity is conserved.

# --- QW-526: SUPERFLUIDITY ---
# Flow past obstacle. Measure drag.

# --- QW-527: NEMATIC ORDER ---
# 2D Grid of spins/phasors.
# Interaction: E = -Sum K(d) cos(theta_i - theta_j)
# Check domain growth.

# --- QW-528: MELTING ---
# Monte Carlo simulation with Temperature T.
# Measure Order Parameter vs T.

# --- QW-529: ELASTICITY ---
# Calculate energy of Twist configuration.
```
