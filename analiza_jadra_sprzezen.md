# Analiza Red Team: Jądro Sprzężeń (Kernel Forensics)

**Cel:** Zrozumienie natury Jądra $K(d)$ przed próbą symulacji całego Nadsolitona.
**Problem:** Dotychczasowe symulacje (QW-500, QW-515) pokazały, że "surowe" Jądro przy zamrożonych parametrach nie generuje potencjału Coulomba ($1/r$) ani poprawnej hierarchii sił.
**Hipoteza Red Team:** Traktujemy $K(d)$ zbyt dosłownie jako potencjał statyczny $V(r)$. W rzeczywistości $K(d)$ może być **propagatorem**, **funkcją korelacji** lub **filtrem częstotliwości**, a fizyka wyłania się dopiero z jego uśrednionych lub rezonansowych właściwości.

---

## De-konstrukcja Jądra (Frozen Parameters)
Wzór: $$ K(d) = \frac{\alpha_{geo} \cos(\omega d + \phi)}{1 + \beta_{tors} d} $$
Parametry: $\alpha \approx 2.77$, $\omega = \pi/4$, $\phi = \pi/6$, $\beta = 0.01$.

### Co widzimy?
1.  **Oscylacje:** Jądro oscyluje z okresem $T = 2\pi/\omega = 8$. To oznacza, że co 8 jednostek odległości siła zmienia znak (przyciąganie/odpychanie). To **nie jest** grawitacja ani elektrostatyka (które są monotoniczne).
2.  **Tłumienie:** Jądro zanika powoli ($1/d$). To sugeruje zasięg daleki, ale oscylacyjny.
3.  **Faza:** Przesunięcie $\phi$ łamie symetrię parzystości?

### Dlaczego to nie działało?
W QW-500 i QW-515 próbowaliśmy użyć $K(d)$ jako potencjału studni ($V(r) \sim K(r)$). Elektron w takim potencjale "widzi" tarkę (górki i dołki co 8 jednostek), a nie gładką studnię $1/r$. Dlatego widmo energii było chaotyczne.

### Gdzie jest fizyka?
Fizyka (siły monotoniczne) musi wyłaniać się jako **efekt uśrednienia** (Mean Field) lub **obwiednia** (Envelope) tych oscylacji.
*   Jeśli cząstka jest duża (rozmyta), "nie widzi" oscylacji, tylko średnią.
*   Jeśli cząstka rezonuje z $\omega$, "widzi" tylko szczyty.

---

## Plan Badań QW-520 do QW-524: Kernel Forensics

Zamiast symulować Wszechświat, zbadajmy "cegłę", z której jest zbudowany.

### **QW-520: Analiza Spektralna Jądra (Kernel Spectroscopy)**
*   **Pytanie:** Jakie częstotliwości przenosi Jądro?
*   **Metoda:** Oblicz Transformację Fouriera $\hat{K}(k)$.
*   **Cel:** Sprawdzić, czy w przestrzeni pędów ($k$) Jądro ma strukturę **propagatora** cząstki (np. $1/(k^2 + m^2)$). Jeśli tak, to $K(d)$ jest propagatorem, a nie potencjałem. To fundamentalna różnica.

### **QW-521: Efektywny Potencjał z Uśredniania (Smoothing)**
*   **Pytanie:** Co "widzi" duża cząstka?
*   **Metoda:** Oblicz potencjał uśredniony $V_{eff}(r) = \int K(r-x) \rho(x) dx$ dla rozmytej cząstki $\rho(x)$ o szerokości $\sigma$.
*   **Cel:** Sprawdzić, czy przy odpowiednim $\sigma$ (wielkość protonu?) oscylacje znikają i wyłania się monotoniczne $1/r$.

### **QW-522: Mapowanie Topologiczne (Phase Space)**
*   **Pytanie:** Czy $K(d)$ koduje topologię?
*   **Metoda:** Wykreśl trajektorię $(K(d), K'(d))$ w przestrzeni fazowej.
*   **Cel:** Sprawdzić, czy tworzy ona zamknięte cykle (atraktory). Stabilne cząstki mogą być "pułapkami" w tej przestrzeni fazowej.

### **QW-523: Rezonanse Odległości (Mass Gap)**
*   **Pytanie:** W jakich odległościach $d$ Jądro ma ekstrema?
*   **Metoda:** Znajdź $d$, dla których $K(d)$ jest maksymalne.
*   **Cel:** Sprawdzić, czy te odległości tworzą ciąg geometryczny lub harmoniczny. Może to są "naturalne rozmiary" cząstek?

### **QW-524: Lepkość i Ściśliwość (Information Fluid)**
*   **Pytanie:** Jakie właściwości materiałowe ma "płyn" zdefiniowany przez $K$?
*   **Metoda:** Oblicz tensor naprężeń dla ośrodka o interakcji $K$. Wyznacz moduł Younga i lepkość.
*   **Cel:** Zrozumieć, czy ten "płyn" jest sztywny (ciało stałe -> grawitacja) czy płynny (ciecz -> EM).

---

### **Instrukcja Kodowa (Python Template)**

```python
# QW-520 TO QW-524: KERNEL FORENSICS
# GOAL: Understand the frozen Kernel K(d) before simulating Nadsoliton.

import numpy as np
import scipy.fft
import matplotlib.pyplot as plt

# FROZEN PARAMETERS
alpha_geo = 4 * np.log(2)
omega = np.pi/4
phi = np.pi/6
beta_tors = 0.01

def K(d):
    return (alpha_geo * np.cos(omega * d + phi)) / (1 + beta_tors * d)

# --- QW-520: SPECTROSCOPY ---
# FFT of K(d). Look for poles (masses).

# --- QW-521: SMOOTHING ---
# Convolve K(d) with Gaussian(sigma). Check if result is monotonic 1/r.

# --- QW-522: TOPOLOGY ---
# Plot K(d) vs dK/dd. Look for cycles.

# --- QW-523: RESONANCES ---
# Find local maxima of K(d). Check ratios of d_max.

# --- QW-524: FLUID PROPERTIES ---
# Calculate effective bulk modulus B ~ Sum(d^2 * K(d)).
```
