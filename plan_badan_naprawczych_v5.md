# Plan Badań QW-540 do QW-544: Evolving Neural Universe

**Paradygmat:** Wszechświat to Ewoluująca Sieć Neuronowa (H9). Prawa fizyki to stan "nauczony" (H11).
**Cel:** Sprawdzić, czy mechanizmy uczenia się (Hebbian Learning) odtwarzają fizykę, której brakowało w modelu "Szkła" (Grawitacja, Czas, Stałe).

---

## Zdefiniowane Zadania (Phase XXIII)

### **QW-540: Grawitacja Hebbowska (Hebbian Gravity)**
*   **Hipoteza:** "Neurons that fire together, wire together". Masa to aktywność.
*   **Cel:** Czy wzmacnianie połączeń między aktywnymi węzłami skraca dystans efektywny ($d_{eff} \propto 1/K$) w sposób przypominający zakrzywienie czasoprzestrzeni?
*   **Metoda:**
    1.  Umieść dwie "masy" (źródła aktywności) w sieci.
    2.  Pozwól połączeniom $K_{ij}$ ewoluować: $\Delta K_{ij} \propto \Psi_i \Psi_j$ (Reguła Hebba).
    3.  Zmierz siłę przyciągania (gradient energii) między masami. Czy jest $1/r^2$?

### **QW-541: Ewolucyjne Dostrojenie (Fine Tuning)**
*   **Hipoteza:** Stałe fizyczne ($\alpha_{geo}, \beta_{tors}$) nie są losowe, ale "nauczone" dla maksymalnego rezonansu.
*   **Cel:** Czy system startujący z losowych parametrów zbiegnie do wartości znanych z FIN ($\alpha \approx 2.77, \beta \approx 0.01$)?
*   **Metoda:**
    1.  Zainicjuj populację Wszechświatów z różnymi parametrami.
    2.  Funkcja celu (Fitness) = Całkowity Rezonans (Energia wiązania).
    3.  Ewolucja (Mutacja + Selekcja).

### **QW-542: Wielki Wybuch i Czas (Arrow of Time)**
*   **Hipoteza:** Czas to proces porządkowania (uczenia się) lub zapominania.
*   **Cel:** Zmierzyć Entropię sieci w czasie ewolucji.
*   **Metoda:**
    1.  Start z chaosu (wysoka entropia).
    2.  Ewolucja Hebbowska.
    3.  Czy Entropia maleje (tworzenie struktur) czy rośnie (zapominanie)? Jak to się ma do II Zasady Termodynamiki?

### **QW-543: Ciemna Energia (Neural Forgetting)**
*   **Hipoteza:** Nieużywane połączenia zanikają (Forgetting). To rozszerza przestrzeń.
*   **Cel:** Czy "zapominanie" w próżni daje efekt odpychania (Ciemna Energia)?
*   **Metoda:**
    1.  Wprowadź regułę zaniku: $\Delta K_{ij} = -\gamma K_{ij}$ (jeśli brak aktywności).
    2.  Zmierz zmianę odległości między spoczynkowymi węzłami.

### **QW-544: Cząstki jako Wspomnienia (Stable Memories)**
*   **Hipoteza:** Cząstki to stabilne atraktory (wspomnienia) w sieci.
*   **Cel:** Czy sieć "pamięta" kształt cząstki po usunięciu wymuszenia?
*   **Metoda:**
    1.  Wymuś wzorzec (np. wir) przez czas $T$.
    2.  Usuń wymuszenie.
    3.  Czy wzorzec przetrwa jako stan własny sieci?

---

### **Instrukcja Kodowa (Python Template)**

```python
# QW-540 TO QW-544: EVOLVING NEURAL UNIVERSE
# PARADIGM: Physics is Learning.

import numpy as np
import matplotlib.pyplot as plt

# --- QW-540: HEBBIAN GRAVITY ---
# Evolve K matrix based on activity.

# --- QW-541: FINE TUNING ---
# Genetic Algorithm for alpha/beta.

# --- QW-542: ARROW OF TIME ---
# Entropy calculation during learning.

# --- QW-543: DARK ENERGY ---
# Decay of connections.

# --- QW-544: PARTICLE MEMORY ---
# Attractor stability check.
```
