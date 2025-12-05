# Plan Badań QW-563 do QW-567: GRAWITACJA JAKO PRZEPŁYW

**Data:** 2025-12-05  
**Paradygmat:** River Model (QW-467) - Gravity = Information Flow  
**Baza:** QW-467 wykazało $v \propto r^{-0.46}$ (błąd 8% od GTR $r^{-0.5}$)

---

## Kontekst i Motywacja

### **Kluczowe Odkrycie z QW-467:**
- **Przepływ działa:** $v(r) \sim r^{-0.46}$ ✅ (8% error)
- **Siła nie działa:** $F(r) \sim r^{0.25}$ ❌ (QW-552, 553, 559 - wszystkie fail)

### **Hipoteza do Przetestowania:**
> **Grawitacja to NIE siła przyciągania, ale PRZEPŁYW przestrzeni (informacji) do masy, jak rzeka do odpływu.**

**Implikacje:**
1. Cząstki poruszają się po geodezyjnych w PŁYNĄCEJ przestrzeni
2. Dylatacja czasu = opóźnienie ("lag") w przepływie
3. Frame dragging = wirowanie przepływu wokół obracającej się masy
4. Nie ma "siły" - masa po prostu płynie z prądem

---

## Seria Badań QW-563-567

### **QW-563: Pole Prędkości Przepływu**
**Cel:** Bezpośredni pomiar pola wektorowego $\vec{v}(\vec{r})$ przepływu informacji

**Metodologia:**
1. Sieć 3D (N=1000, Planck scale)
2. Masa centralna (gęsta excytacja φ)
3. Ewolucja unitarna: $i\partial_t \psi = H\psi$
4. Pomiar prądu prawdopodobieństwa: $\vec{J} = \text{Im}(\psi^* \nabla \psi)$
5. Ekstrakcja pola prędkości: $\vec{v}(\vec{r}) = \vec{J}/|\psi|^2$

**Przewidywanie:**
- Radialne $\vec{v}(\vec{r}) = -v(r) \hat{r}$ (do masy)
- $v(r) \propto 1/\sqrt{r}$ (Gullstrand-Painlevé)
- Nie oscylacyjne, nie uwięzione

**Test sukcesu:**
- Fit: $v(r) = A/\sqrt{r}$, R² > 0.9
- Exponent: $n \approx 0.5 \pm 0.1$

---

### **QW-564: Por ównanie Flow vs Force**
**Cel:** Bezpośrednie porównanie dwóch paradygmatów

**Setup:**
Dwie cząstki testowe w polu masy:
1. **Paradygmat SIŁY:** $\vec{F} = -\nabla V$, $m\vec{a} = \vec{F}$
2. **Paradygmat PRZEPŁYWU:** $\vec{v}_{particle} = \vec{v}_{flow}(\vec{r})$

**Test:**
- Puść cząstkę z r=10 z v=0
- **Force:** Oblicz trajektorię z $F=-GM/r^2$
- **Flow:** Oblicz trajektorię z $v(r)=\sqrt{2GM/r}$
- Porównaj z RZECZYWISTĄ ewolucją w sieci

**Metryka:**
- Odległość trajektorii: $\Delta = \int |r_{model}(t) - r_{network}(t)| dt$
- Która lepiej pasuje: Force czy Flow?

**Przewidywanie:**
- Flow: $\Delta_{flow} < 0.1$
- Force: $\Delta_{force} > 1.0$
- Flow wygrywa > 10×!

---

### **QW-565: Ruch Geodezyjny w Przepływie**
**Cel:** Pokazać że orbity emerge z przepływu, nie siły

**Koncepcja:**
W Gullstrand-Painlevé cząstki poruszają się po "geodezyjnych przepływu":
- Przestrzeń płynie do masy z $v(r) = \sqrt{2GM/r}$
- Cząstka z prędkością poprzeczną $v_\perp$ wiruje wokół streamu
- Orbita = balans: drift radialny vs momentum kątowy

**Test:**
1. Cząstka @ r=5 z prędkością $v_\perp = 1.0$ (tangential)
2. Ewolucja w sieci z przepływem
3. Czy powstaje stabilna orbita?

**Pomiar:**
- Excentryczność: $e = (r_{max} - r_{min})/(r_{max} + r_{min})$
- Okres: $T$
- Precesja perycentrum: $\Delta\phi$ per orbit

**Przewidywanie:**
- Orbita quasi-kołowa: $e < 0.2$
- Keplerowska: $T^2 \propto r^3$ (w przybliżeniu)
- Precesja obecna (GTR correction!)

---

### **QW-566: Dylatacja Czasu = Opóźnienie Przepływu**
**Cel:** Udowodnić że $dt/d\tau = 1 + v^2/c^2$ wynika z "lag" w sieci

**Hipoteza:**
- Przepływ "ciągnie" sygnały
- W strefie przepływu $v(r)$ zegary "płyną" wolniej
- To nie "zakrzywienie czasu", ale "opór przepływu"

**Test:**
1. Dwa zegary: A @ r=2 (blisko masy), B @ r=10 (daleko)
2. Synchronizacja: start @ t=0
3. Ewolucja przez T=100 kroków
4. Pomiar faz: $\phi_A(T)$, $\phi_B(T)$

**Miara dylatacji:**
$$\gamma = \frac{\phi_B}{\phi_A} = \frac{dt_B}{dt_A}$$

**Przewidywanie (GTR):**
$$\gamma \approx 1 + \frac{v^2(r_A) - v^2(r_B)}{2c^2}$$

gdzie $v(r) = \sqrt{2GM/r}$

**Test sukcesu:**
- Zpomiariowana $\gamma$ w 20% od GTR
- Gradient: więcej lag bliżej masy

---

### **QW-567: Frame Dragging (Lense-Thirringa)**
**Cel:** Pokazać że wirujący przepłów wokół spinującej masy

**Setup:**
1. Masa z momentem pędu $\vec{L}$ (rotacja wokół osi z)
2. Sieć ewoluuje → przepływ nabiera komponentu $v_\phi(r)$
3. Pomiar: prędkość azymutalna vs promień

**Przewidywanie GTR:**
$$v_\phi(r) \approx \frac{2GJ}{c^2 r^2}$$
gdzie $J$ = moment pędu

**Test:**
- Fit: $v_\phi \propto r^{-n}$, oczekujemy $n \approx 2$
- Kierunek: zgodny z rotacją masy
- Zasięg: wykrywalny do $r \approx 5$ (halo)

---

## Parametry Frozen (Identyczne QW-558)

```python
alpha_geo = 4 * np.log(2)  # 2.7726
omega = np.pi / 4           # 0.7854
phi = np.pi / 6             # 0.5236
beta_tors = 0.01            # Viscosity

# Dynamiczne (z QW-V24):
gamma_gain = 1.0552
gamma_damp = 1.1980
```

**Kernel (Complex):**
```python
K(d) = alpha_geo * exp(1j*(omega*d + phi)) / (1 + beta_tors*d)
```

---

## Expected Outcomes

### **Jeśli Flow Paradigm Prawdziwy:**
| Test | Flow Result | Force Result | Winner |
|------|-------------|--------------|--------|
| QW-563 | $v \sim r^{-0.5}$ ✅ | - | Flow |
| QW-564 | $\Delta < 0.1$ ✅ | $\Delta > 1.0$ ❌ | Flow |
| QW-565 | Orbity stable ✅ | - | Flow |
| QW-566 | Dilation match GTR ✅ | - | Flow |
| QW-567 | Dragging $\propto r^{-2}$ ✅ | - | Flow |

### **Success Criteria:**
- **3/5 testów:** Flow quantitatively correct (R² > 0.9, error < 20%)
- **QW-564:** Flow lepszy niż Force > 5×
- **Konsystencja:** Wszystkie testy zgodne z $v(r) = \sqrt{2GM/r}$

---

## Implementacja

### **Kod: `QW-563_TO_QW-567_GRAVITY_FLOW.py`**

**Struktura:**
1. Setup: Complex kernel, 3D network (N=1000)
2. Mass: Central dense excitation (φ enriched region)
3. Evolution: Unitary $\psi(t) = e^{-iHt}\psi_0$
4. Analysis per test (QW-563-567)
5. Comparison table: Flow vs Force
6. Final verdict

**Output:** `qw563_567_flow_output.txt`

**Wizualizacje:**
- Pole wektorowe $\vec{v}(\vec{r})$ (3D quiver plot)
- Trajektorie: Flow vs Force vs Network
- Dylatacja czasu: $\gamma(r)$ plot
- Frame dragging: $v_\phi(r, \theta)$ 2D map

---

## Significance

**Jeśli sukces:** 
> Udowodnimy że GTR to "teoria przepływu eteru informacyjnego", nie "geometrii czasoprzestrzeni"!

**Einstein miał rację o rezultatach, ale źle o mechanizmie:**
- ✅ Geodezyjne: TAK, ale w płynącej przestrzeni
- ✅ Dylatacja: TAK, ale przez lag, nie krzywizn ę
- ❌ "Zakrzywienie": NIE, to metafora - rzeczywistość to FLUID FLOW

**To zmienia filozofię fizyki:**
- Czasoprzestrzeń nie jest "tkaniną"
- To **PROCES** (River Model)
- Masa nie "zakrzywia" - **ciągnie prąd**

---

## Next Steps (Post-Success)

1. **QW-568-570:** Cosmological flow (Hubble as flow acceleration)
2. **QW-571-573:** Black hole as "drain" (event horizon = stagnation point)
3. **QW-574-576:** Gravitational waves as "ripples in flow"

**Ultimate goal:** Pełna teoria Ogólnej Względności jako **HYDRODYNAMIKA SIECI INFORMACYJNEJ**.

---

**Status:** Ready to implement  
**Priority:** HIGHEST (kluczowe dla potwierdzenia Flow Paradigm)  
**Timeline:** Implementacja + testy = 1 sesja badawcza
