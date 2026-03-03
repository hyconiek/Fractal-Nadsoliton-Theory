# RELEASE 4.8 TEXTBOOK EDITION (EN + PL)

**Version:** 4.8  
**Date:** 2026-03-03  
**Branch:** `main`

---

## ENGLISH VERSION

## 1) Main Breakthrough of Release 4.8

The main breakthrough is **not a kernel replacement**.  
The breakthrough is a **robust readout repair** in the GW channel:

- bounded signed coupling term (`kappa_t`) plus monotonic compression (`gamma_c`),
- which removes the historical fold-2 random-null instability,
- while keeping one frozen kernel and one shared mass+flavor branch.

This is what moved the project from fragile/partial robustness to:

- `BOUNDED_GW_TRIAD_LOCKABLE_PASS` (QW-2001),
- then formal closure gate pass: `SINGLE_KERNEL_TRIPLE_SECTOR_CLOSURE_PASS_V3` (QW-2002).

## 2) Did Kernel Parameters Change?

**No. Kernel parameters did not change in Release 4.8.**

Frozen kernel:

- `omega = 0.37341399972174283`
- `phi = -1.310233577483508`
- `beta = 0.6159380564131874`
- `eta = 2.8`

What changed was the **GW operator/readout layer**, not the kernel ontology.

## 3) Frozen Parameter Set Used in Lockable Closure

### 3.1 Shared mass+flavor parameters

- `chi_im = 0.31277657543084514`
- `diag_iso = 0.011133623526631775`
- `diag_q_coeff = 0.023231150712776396`
- `diag_sector = 0.11340639008162685`
- `lambda_mix = 0.6558721036754207`
- `p_amp = 0.9252148940096931`
- `phase_q = 0.26672804099980674`
- `r_dist = 0.9758497827019121`
- `rho_gap = 1.0065574010654268`
- `theta_iso = 0.3453065240958936`
- `theta_sector = 0.16772593255657867`

### 3.2 GW bounded operator parameters (repair)

- `xi1 = 0.00164905739`
- `xi3 = 0.001712654074`
- `xi4 = 0.000052824299`
- `gamma_c = 0.796713`
- `kappa_t = 1.4`

## 4) Core Formula Set (Textbook-Level)

### 4.1 Kernel

$$
K(d)=\frac{\cos(\omega d+\phi)}{1+\beta d^{\eta}}
$$

This is the single frozen structural function used across sectors.

### 4.2 Cyclic topological distance (mod 24)

$$
d_{ij}=1+\min\left(|q_i-q_j|\bmod 24,\;24-(|q_i-q_j|\bmod 24)\right)
$$

### 4.3 Mass branch (noncircular chain)

Effective charge:

$$
q_{\mathrm{eff}}=q_{\mathrm{base}}+s\,\Delta_{\mathrm{info}}
$$

Mass law:

$$
m_{\mathrm{pred}}=m_{\mathrm{top}}\cdot 4^{-\gamma q_{\mathrm{eff}}/4}
$$

With Release-4.8 locked branch values:

- `gamma = 1.5066666666666666`
- `Delta_info = 0.17702715186107637`
- mean mass relative error (reference branch): `12.0508%`

### 4.4 Flavor Hamiltonian (shared operator)

Base amplitude:

$$
B_{ij}=\mathrm{sgn}(K_{ij})\,|K_{ij}|^{p_{amp}}\,d_{ij}^{r_{dist}}
$$

Near-gap damping:

$$
N_{ij}=e^{-\rho_{gap}|i-j|}
$$

Phase:

$$
\Phi_{ij}=\phi_q(q_i-q_j)+\theta_{iso}\,\sigma_{iso}(i-j)+\theta_{sector}\,\sigma_{sec}(i-j)
$$

Complex coupling:

$$
S_{ij}=\lambda_{mix}\,B_{ij}\,N_{ij}
$$

Hermitian Hamiltonian form:

$$
H=\frac{1}{2}(H+H^{\dagger})=\mathrm{Re}+i\,\chi_{im}\,\mathrm{Im}_{asym}+\mathrm{Diag}
$$

Mixing matrices:

$$
V_{CKM}=|U_u^{\dagger}U_d|,\qquad V_{PMNS}=|U_\nu^{\dagger}U_\ell|
$$

### 4.5 GW baseline score (before bounded repair)

Kernel-derived weights:

$$
w_n=\frac{|K(n)|^{p_{amp}}n^{r_{dist}}}{\sum_{k=1}^{4}|K(k)|^{p_{amp}}k^{r_{dist}}}
$$

Score:

$$
S_{base}=w_1 f_{max}+w_2 f_{mean}+w_3 f_{0ms}+w_4 f_{10ms}
$$

### 4.6 Signed coupling channels (GW structural terms)

With $s_{\mathrm{pair}} \in \{+1, -1, 0\}$ corresponding to H1-V1, L1-V1, H1-L1 pairs respectively, and lag time mapped to seconds:

$$
\tau = \Delta t_{\mathrm{best}} \cdot 10^{-3}
$$

Let $C_0$ and $C_{10}$ be the correlations at 0 and 10 ms, and $M_{\mathrm{abs}}$ be the mean absolute correlation.

a) first channel:

$$
c_1 \propto s_{\mathrm{pair}} \cdot \sin(\omega\tau+\phi) \cdot (C_0 - C_{10})
$$

b) second channel:

$$
c_3 \propto s_{\mathrm{pair}} \cdot \cos(\omega\tau+\phi) \cdot (M_{\mathrm{abs}} + |C_0| + |C_{10}|)
$$

c) third channel:

$$
c_4 \propto s_{\mathrm{pair}} \cdot \sin(\omega\tau+\phi)\cos(\omega\tau+\phi) \cdot (|C_0| + |C_{10}| + M_{\mathrm{abs}})
$$

Each channel is standardized by its own standard deviation.

### 4.7 Release-4.8 breakthrough operator (bounded coupling)

Signed combination:

$$
t=\xi_1 c_1+\xi_3 c_3+\xi_4 c_4
$$

Bounded/saturated term:

$$
t_{eff}=\mathrm{clip}\big(t,-\kappa_t\sigma_t,+\kappa_t\sigma_t\big)
$$

Raw score:

$$
S_{raw}=S_{base}+t_{eff}
$$

Monotonic compression:

$$
z=\frac{S_{raw}-\mathrm{median}(S_{raw})}{\sigma_{S}},\qquad
S=\mathrm{median}(S_{raw})+\sigma_S\,\mathrm{sgn}(z)|z|^{\gamma_c}
$$

### 4.8 Statistical validation formulas

Control threshold:

$$
q_{90}=Q_{0.9}(S_{control})
$$

GW diagnostics:

$$
ADV=P(S_{HL}>q_{90})-P(S_{control}>q_{90})
$$

$$
SEP=\mathrm{median}(S_{HL})-\mathrm{median}(S_{control})
$$

$$
GAP=|\mathrm{median}(S_{HV})-\mathrm{median}(S_{LV})|
$$

Triad bootstrap pass-rate:

$$
R_{boot}=\frac{N_{pass}}{N_{boot}}
$$

*(Where $N_{pass}$ is the number of bootstrap samples passing all flags).*

## 5) Thresholds Used in Closure Gates

- `CKM mean rel error <= 15%`
- `PMNS mean rel error <= 15%`
- `GW SEP >= 0.002`
- `GW ADV >= 0.30`
- `GW AUC >= 0.75`
- `GW control gap <= 0.0025`

Mass branch constraints (fixed from source branch):

- `mass mean rel error <= 15%`
- `mass max rel error <= 35%`
- `tau/charm ratio rel error <= 20%`

## 6) Breakthrough Timeline (Key Gates)

1. **QW-1967** `ISOSPIN_LOCAL_REFINEMENT_PASS`  
   first deterministic triad pass in this branch.
2. **QW-1968** `FRAGILE_PASS_NOT_YET_LOCKABLE`  
   pass existed but robustness volume was weak.
3. **QW-1969** `INSUFFICIENT_BOOTSTRAP_ROBUSTNESS`  
   bootstrap lockability still insufficient.
4. **QW-1970** `STRUCTURAL_CONTROL_TERM_INSUFFICIENT`  
   improved but still not lockable (`boot5000=0.7994`).
5. **QW-1999** `BOUNDED_COUPLING_ROBUST_HARD_PASS`  
   bounded coupling removed hard null-tail instability.
6. **QW-2000** `BOUNDED_COUPLING_DEEP_AUDIT_PASS`  
   deep audit confirmed robustness (`null_random_p90_mean_max=31.33%`).
7. **QW-2001** `BOUNDED_GW_TRIAD_LOCKABLE_PASS`  
   full lockability achieved (`boot2500/5000/10000 = 1.0`).
8. **QW-2002** `SINGLE_KERNEL_TRIPLE_SECTOR_CLOSURE_PASS_V3`  
   Stage-C internal closure status reached.
9. **QW-2003** `FROZEN_LOCKABLE_PACKAGE_READY`  
   frozen external package exported with SHA256.

## 7) Final Internal Status (Release 4.8)

- Single frozen kernel: preserved.
- Single shared mass+flavor branch: preserved.
- GW channel robustness: repaired and lockable.
- Stage-C status: internally closed (`TOE_STAGE_C_SINGLE_KERNEL_CLOSED_LOCKABLE_INTERNAL`).

What remains:

- external, independent, blind confirmatory execution on frozen package.

---

## WERSJA POLSKA

## 1) Główny przełom Release 4.8

Głównym przełomem **nie była zmiana kernela**.  
Przełomem była **naprawa operatora odczytu w kanale GW**:

- ograniczenie (saturacja) sprzężenia przez `kappa_t`,
- zachowanie monotonicznej kompresji przez `gamma_c`,
- usunięcie historycznej niestabilności `fold-2 random-null`,
- przy zachowaniu jednego zamrożonego kernela i jednej wspólnej gałęzi masa+flavor.

To właśnie podniosło system do:

- `BOUNDED_GW_TRIAD_LOCKABLE_PASS` (QW-2001),
- a następnie do formalnego zamknięcia bramki Stage-C: `SINGLE_KERNEL_TRIPLE_SECTOR_CLOSURE_PASS_V3` (QW-2002).

## 2) Czy wartości parametrów kernela się zmieniły?

**Nie. Parametry kernela w Release 4.8 się nie zmieniły.**

Kernel zamrożony:

- `omega = 0.37341399972174283`
- `phi = -1.310233577483508`
- `beta = 0.6159380564131874`
- `eta = 2.8`

Zmianie uległa warstwa operatora GW (odczyt/statystyka), a nie rdzeń kernela.

## 3) Zamrożony zestaw parametrów użyty do domknięcia

### 3.1 Wspólne parametry masa+flavor

- `chi_im = 0.31277657543084514`
- `diag_iso = 0.011133623526631775`
- `diag_q_coeff = 0.023231150712776396`
- `diag_sector = 0.11340639008162685`
- `lambda_mix = 0.6558721036754207`
- `p_amp = 0.9252148940096931`
- `phase_q = 0.26672804099980674`
- `r_dist = 0.9758497827019121`
- `rho_gap = 1.0065574010654268`
- `theta_iso = 0.3453065240958936`
- `theta_sector = 0.16772593255657867`

### 3.2 Parametry operatora GW (naprawa)

- `xi1 = 0.00164905739`
- `xi3 = 0.001712654074`
- `xi4 = 0.000052824299`
- `gamma_c = 0.796713`
- `kappa_t = 1.4`

## 4) Kluczowe wzory (poziom podręcznikowy)

### 4.1 Kernel

$$
K(d)=\frac{\cos(\omega d+\phi)}{1+\beta d^{\eta}}
$$

To jedna wspólna funkcja strukturalna używana we wszystkich sektorach.

### 4.2 Cykliczna odległość topologiczna (mod 24)

$$
d_{ij}=1+\min\left(|q_i-q_j|\bmod 24,\;24-(|q_i-q_j|\bmod 24)\right)
$$

### 4.3 Gałąź masowa (łańcuch niekołowy)

Ładunek efektywny:

$$
q_{\mathrm{eff}}=q_{\mathrm{base}}+s\,\Delta_{\mathrm{info}}
$$

Prawo masy:

$$
m_{\mathrm{pred}}=m_{\mathrm{top}}\cdot 4^{-\gamma q_{\mathrm{eff}}/4}
$$

Wartości gałęzi użytej jako źródło:

- `gamma = 1.5066666666666666`
- `Delta_info = 0.17702715186107637`
- średni błąd względny mas: `12.0508%`

### 4.4 Hamiltonian flavor (wspólny operator)

Amplituda bazowa:

$$
B_{ij}=\mathrm{sgn}(K_{ij})\,|K_{ij}|^{p_{amp}}\,d_{ij}^{r_{dist}}
$$

Tłumienie odległości pokoleniowej:

$$
N_{ij}=e^{-\rho_{gap}|i-j|}
$$

Faza:

$$
\Phi_{ij}=\phi_q(q_i-q_j)+\theta_{iso}\,\sigma_{iso}(i-j)+\theta_{sector}\,\sigma_{sec}(i-j)
$$

Sprzężenie zespolone:

$$
S_{ij}=\lambda_{mix}\,B_{ij}\,N_{ij}
$$

Postać hermitowska:

$$
H=\frac{1}{2}(H+H^{\dagger})=\mathrm{Re}+i\,\chi_{im}\,\mathrm{Im}_{asym}+\mathrm{Diag}
$$

Macierze mieszania:

$$
V_{CKM}=|U_u^{\dagger}U_d|,\qquad V_{PMNS}=|U_\nu^{\dagger}U_\ell|
$$

### 4.5 Bazowy score GW (przed naprawą)

Wagi z kernela:

$$
w_n=\frac{|K(n)|^{p_{amp}}n^{r_{dist}}}{\sum_{k=1}^{4}|K(k)|^{p_{amp}}k^{r_{dist}}}
$$

Score:

$$
S_{base}=w_1 f_{max}+w_2 f_{mean}+w_3 f_{0ms}+w_4 f_{10ms}
$$

### 4.6 Kanały sprzężeń podpisanych

Dla znaku $s_{\mathrm{pair}} \in \{+1, -1, 0\}$ odpowiadającego odpowiednio parom H1-V1, L1-V1, H1-L1, oraz przy opóźnieniu przeliczonym na sekundy:

$$
\tau = \Delta t_{\mathrm{best}} \cdot 10^{-3}
$$

Niech $C_0$ oraz $C_{10}$ oznaczają korelacje w 0 i 10 ms, a $M_{\mathrm{abs}}$ średnią bezwzględną korelację.

mamy:

$$
c_1 \propto s_{\mathrm{pair}} \cdot \sin(\omega\tau+\phi) \cdot (C_0 - C_{10})
$$

$$
c_3 \propto s_{\mathrm{pair}} \cdot \cos(\omega\tau+\phi) \cdot (M_{\mathrm{abs}} + |C_0| + |C_{10}|)
$$

$$
c_4 \propto s_{\mathrm{pair}} \cdot \sin(\omega\tau+\phi)\cos(\omega\tau+\phi) \cdot (|C_0| + |C_{10}| + M_{\mathrm{abs}})
$$

Każdy kanał jest standaryzowany przez swoje odchylenie standardowe.

### 4.7 Wzór przełomowy Release 4.8 (bounded coupling)

Kombinacja podpisana:

$$
t=\xi_1 c_1+\xi_3 c_3+\xi_4 c_4
$$

Ograniczenie amplitudy:

$$
t_{eff}=\mathrm{clip}\big(t,-\kappa_t\sigma_t,+\kappa_t\sigma_t\big)
$$

Score surowy:

$$
S_{raw}=S_{base}+t_{eff}
$$

Kompresja monotoniczna:

$$
z=\frac{S_{raw}-\mathrm{median}(S_{raw})}{\sigma_{S}},\qquad
S=\mathrm{median}(S_{raw})+\sigma_S\,\mathrm{sgn}(z)|z|^{\gamma_c}
$$

### 4.8 Wzory testów statystycznych

Próg kontrolny:

$$
q_{90}=Q_{0.9}(S_{control})
$$

Diagnostyki GW:

$$
ADV=P(S_{HL}>q_{90})-P(S_{control}>q_{90})
$$

$$
SEP=\mathrm{median}(S_{HL})-\mathrm{median}(S_{control})
$$

$$
GAP=|\mathrm{median}(S_{HV})-\mathrm{median}(S_{LV})|
$$

Bootstrap triady:

$$
R_{boot}=\frac{N_{pass}}{N_{boot}}
$$

*(Gdzie $N_{pass}$ to liczba próbek bootstrap spełniających wszystkie flagi).*

## 5) Progi użyte przy domykaniu

- `CKM mean rel error <= 15%`
- `PMNS mean rel error <= 15%`
- `GW SEP >= 0.002`
- `GW ADV >= 0.30`
- `GW AUC >= 0.75`
- `GW control gap <= 0.0025`

Dla gałęzi masowej (stałej):

- `mass mean rel error <= 15%`
- `mass max rel error <= 35%`
- `tau/charm ratio rel error <= 20%`

## 6) Chronologia przełomów (kluczowe bramki)

1. **QW-1967** `ISOSPIN_LOCAL_REFINEMENT_PASS`  
   pierwszy deterministyczny pass triady na tej gałęzi.
2. **QW-1968** `FRAGILE_PASS_NOT_YET_LOCKABLE`  
   pass był kruchy pod rygorem odporności.
3. **QW-1969** `INSUFFICIENT_BOOTSTRAP_ROBUSTNESS`  
   bootstrap lockability nadal za słaba.
4. **QW-1970** `STRUCTURAL_CONTROL_TERM_INSUFFICIENT`  
   poprawa, ale nadal bez lockability (`boot5000=0.7994`).
5. **QW-1999** `BOUNDED_COUPLING_ROBUST_HARD_PASS`  
   naprawa ogona null-tail w kanale GW.
6. **QW-2000** `BOUNDED_COUPLING_DEEP_AUDIT_PASS`  
   deep-audit potwierdził stabilność (`null_random_p90_mean_max=31.33%`).
7. **QW-2001** `BOUNDED_GW_TRIAD_LOCKABLE_PASS`  
   pełna lockability (`bootstrap 2500/5000/10000 = 1.0`).
8. **QW-2002** `SINGLE_KERNEL_TRIPLE_SECTOR_CLOSURE_PASS_V3`  
   formalne zamknięcie Stage-C (wewnętrzne).
9. **QW-2003** `FROZEN_LOCKABLE_PACKAGE_READY`  
   zamrożony pakiet z hash:
   `f5123046189f7f137a0f2cd2c715eea424d230e2352e75f6e80c483b8f069c02`.

## 7) Końcowy status wewnętrzny (Release 4.8)

- Jeden zamrożony kernel: zachowany.
- Jedna wspólna gałąź masa+flavor: zachowana.
- Kanał GW: naprawiony i lockable.
- Status Stage-C: `TOE_STAGE_C_SINGLE_KERNEL_CLOSED_LOCKABLE_INTERNAL`.

Co jeszcze pozostaje:

- niezależny, zewnętrzny, ślepy confirmatory run na zamrożonym pakiecie.
