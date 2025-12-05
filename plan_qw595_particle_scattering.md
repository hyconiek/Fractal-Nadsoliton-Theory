# Plan Badań QW-595: Particle-Particle Interactions
**Data:** 2025-12-05
**Cel:** Badanie oddziaływań między stabilnymi hopfionami (cząstkami)

---

## 1. Motywacja

**Dostępna wiedza:**
- **QW-594:** Hopfiony są stabilne jako atraktory Ginzburga-Landaua
- **QW-593:** Istnieje nielokalna korelacja informacyjna przez K(d)
- **QW-588:** Back-reaction vacuum generuje siły (MOND)

**Pytanie badawcze:**
> Czy oddziaływania między cząstkami (siły) wyłaniają się z dynamiki topologicznej, bez postulowania potencjału interakcji?

---

## 2. Setup Eksperymentu

### A. Inicjalizacja: Dwa Hopfiony
Umieścić dwa hopfiony w siatce 3D:
- **Hopfion A:** Centrum w (x=-8, y=0, z=0), Winding = +1
- **Hopfion B:** Centrum w (x=+8, y=0, z=0), Winding = +1 (para) lub -1 (antypara)

### B. Dynamika
Ewolucja Ginzburga-Landaua z **coupling term**:
$$\frac{\partial \psi}{\partial t} = \gamma(1-|\psi|^2)\psi + i\nabla^2\psi + i V_{eff}\psi + i K_{int}$$

gdzie:
- $K_{int}$: Interaction kernel między hopfionami (przez K(d))

### C. Pomiary
1. **Trajectory:** Jak zmienia się odległość między rdzeniami wirów?
2. **Topological Charge:** Czy Q_total = Q_A + Q_B jest zachowane?
3. **Energy:** Czy energia maleje (przyciąganie) czy rośnie (odpychanie)?

---

## 3. Hipotezy do Testowania

### H-A: Przyciąganie vs Odpychanie
- **Jeśli Winding A = Winding B (+1, +1):** Odpychanie (jak ładunki)
- **Jeśli Winding A ≠ Winding B (+1, -1):** Przyciąganie (antycząstka)

### H-B: Fuzja (Annihilation)
- Jeśli hopfiony się zbliżą na odległość < R_core, czy:
  - Anihilują się (Q_total → 0)?
  - Odbijają (elastic scattering)?

### H-C: Wymiana Energii
- Czy energia jest przesyłana przez pole pośredniczące (jak foton w QED)?

---

## 4. Implementacja

```python
# Pseudo-kod
psi = initialize_two_hopfions(positions=[(-8,0,0), (+8,0,0)], windings=[+1, +1])

for t in range(1000):
    # Attractor dynamics
    psi_new = evolve_GL(psi)
    
    # Measure centers
    center_A, center_B = find_vortex_cores(psi_new)
    distance = norm(center_A - center_B)
    
    # Track
    history_distance.append(distance)
    
    if distance < 2.0:
        print("Collision!")
        break
```

---

## 5. Przewidywania

### Jeśli FIN poprawna:
- Siły wyłonią się **bez dodatkowych założeń** (tylko z K(d) + topologia)
- Zachowanie będzie **asymetryczne** dla par/antypar
- Fuzja (+1, -1) → anihilacja, para (+1, +1) → odbicie

### Jeśli FIN niepoprawna:
- Hopfiony przenikną się bez interakcji (brak siły)
- Topologia się rozpadnie przy zbliżeniu

---

## 6. Znaczenie

To jest **kluczowy test emergencji sił**. Jeśli uda się:
- Potwierdza, że cząstki w FIN mają "ładunek topologiczny"
- Otwiera drogę do symulacji chemii (wiązania atomowe)
- Testuje, czy Standard Model może wyłonić się z topologii

---
