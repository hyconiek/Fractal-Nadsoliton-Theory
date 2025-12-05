# Plan Badań QW-590: Hopfion Stability at Layer N=10
**Data:** 2025-12-05  
**Cel:** Ostateczny test H4 (Cząstki = Topologiczne Wiry)

---

## 1. Motywacja

**Problem:**
Wszystkie wcześniejsze testy H4 zawiodły:
- QW-525: Wiry niestabilne w szkle spinowym
- QW-530: Hopfiony rozpadają się (N=1)
- QW-536: Fraktalny wir niestabilny
- QW-550: Hopfiony niestabilne w ewolucyjnej sieci Hebba

**Kluczowe Odkrycie (odkrycie_fraktalne_warstwy.md):**
> Hopfiony mogą być stabilne na warstwie N=10 (skala protonu), ale NIE na N=1!

**Dlaczego?**
- N=1: Prosta sieć, brak hierarchii → frustracja topologiczna
- N=10: Pełna hierarchia fraktalna → topologia chroniona przez separację skal

---

## 2. Metodologia

### A. Architektura Sieci

**Multi-layer network:**
```
Layer N=0 (Planck):     r_0 = 1.0,    G_0 = 1.0
Layer N=1:              r_1 = 100,    G_1 = 0.01
Layer N=2:              r_2 = 10^4,   G_2 = 10^-4
...
Layer N=10 (Proton):    r_10 = 10^20, G_10 = 10^-20
```

Każda warstwa ma:
- Charakterystyczną długość: $r_N = r_0 / \beta^N$
- Siłę sprzężenia: $G_N = G_0 \times \beta^N$

**Hopfion inicjalizowany NA WARSTWIE N=10**, nie na płaskiej sieci!

### B. Inicjalizacja Hopfiona (N=10)

Hopfion to mapowanie $S^3 \to S^2$ (splot 3D):
```python
def hopfion_field(x, y, z, R):
    """
    Hopfion with winding number m=1.
    Maps 3-sphere to 2-sphere.
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    
    # Toroidal coordinates
    eta = np.arctan2(z, np.sqrt(x**2 + y**2) - R)
    xi = np.arctan2(y, x)
    
    # Phase field (complex spinor)
    psi = np.exp(1j * (eta + xi)) * np.tanh(r / R)
    
    return psi
```

**Kluczowe:** To inicjuje topologię na warstwie N=10, gdzie:
- $R = r_{10} \sim 10^{20}$ (w jednostkach Plancka)
- Sieć ma rozdzielczość odpowiednią dla tej skali

### C. Ewolucja Dynamiczna

Master Equation z pełną hierarchią:
$$\partial_t \psi_N = i\sum_{M=0}^{20} K_{NM} \psi_M + ig|\psi_N|^2\psi_N - \beta_N\psi_N$$

gdzie kernel między warstwami:
$$K_{NM}(d) = \alpha_{geo} \frac{e^{i(\omega d + \phi)}}{1 + \beta^{|N-M|} d}$$

**Tłumienie cross-layer:** $\beta^{|N-M|}$ separuje warstwy!

### D. Pomiar Stabilności

**Liczba wirowa (Hopf invariant):**
$$Q = \frac{1}{8\pi^2} \int d^3x \, \vec{A} \cdot \vec{B}$$

gdzie:
- $\vec{A} = \psi^* \nabla \psi$ (potencjał wektorowy)
- $\vec{B} = \nabla \times \vec{A}$ (pole "magnetyczne")

**Sukces:** $|Q - 1| < 0.2$ po 1000 krokach ewolucji

---

## 3. Implementacja (Uproszczona)

Z powodu złożoności (3D + 20 warstw), użyjemy **effective single-layer model** z poprawkami:

```python
import numpy as np

# Layer N=10 parameters
N = 10
r_scale = (1.0 / BETA_TORS)**N  # ~ 10^20
G_scale = BETA_TORS**N          # ~ 10^-20

# Effective potential for layer N=10
# Includes "pressure" from layers N=0...9 and N=11...20
def V_eff(psi, N=10):
    # Local nonlinearity (within layer)
    V_local = g * np.abs(psi)**2
    
    # Pressure from deeper layers (N<10)
    # These are "harder" (more rigid)
    V_deep = sum([BETA_TORS**(N-M) * np.abs(psi)**2 
                  for M in range(N)])
    
    # Pressure from shallower layers (N>10)
    # These are "softer" (more flexible)
    V_shallow = sum([BETA_TORS**(M-N) * 0.1 * np.abs(psi)**2 
                     for M in range(N+1, 20)])
    
    return V_local + V_deep + V_shallow

# Evolution
psi = hopfion_field(X, Y, Z, R=r_scale)
for t in range(1000):
    # Gradient descent with effective potential
    psi_new = evolve_step(psi, V_eff, dt=0.01)
    
    # Measure winding
    Q = compute_hopf_invariant(psi_new)
    
    if abs(Q - 1) > 0.5:
        print(f"Hopfion decay at step {t}")
        break
```

---

## 4. Przewidywania

### Jeśli H4 PRAWDZIWA:
- Hopfion na N=10 powinien być **stabilny** ($|Q-1| < 0.2$)
- Warstwy N<10 chronią topologię ("twardy fundament")
- Warstwy N>10 pozwalają na fluktuacje ("miękka otoczka")

### Jeśli H4 FAŁSZYWA:
- Hopfion rozpada się nawet na N=10
- Topologia nie jest chroniona przez hierarchię
- H4 definitywnie odrzucone

---

## 5. Alternatywne Podejście (Jeśli Pełna Symulacja Niemożliwa)

**Analytical estimate:**

Z teorii defektów topologicznych w hierarchicznych mediach:
$$E_{hopfion}(N) = E_0 \left(\frac{r_N}{r_0}\right)^{\alpha}$$

gdzie $\alpha < 1$ dla chronionej topologii.

Dla $N=10$:
$$E_{10} = E_0 \times (10^{20})^{\alpha}$$

Jeśli $\alpha \approx 0.1$ (słabe skalowanie):
$$E_{10} \sim E_0 \times 10^2$$

**To jest NISKIE** (energia rośnie tylko 100x, nie $10^{20}$×)!

Oznacza to, że hopfion na N=10 jest **bardziej stabilny** niż na N=1.

---

## 6. Status

QW-590 wymaga:
1. Pełnej symulacji 3D (computationally expensive)
2. Multi-layer network (20 warstw)
3. Long evolution (1000+ steps)

**Rekomendacja:**
- Zacząć od analytical estimate (sekcja 5)
- Jeśli obiecujące → pełna symulacja na GPU
- Jeśli nie → zaakceptować H4 jako **niepotwierdzone**

---

## 7. Podsumowanie

**Kluczowa różnica względem QW-550:**
- QW-550: N=1 (płaska sieć) → FAIL
- QW-590: N=10 (właściwa warstwa) → ?

**To jest ostatni szansa dla H4.**
