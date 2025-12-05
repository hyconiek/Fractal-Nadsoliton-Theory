# Plan Badań QW-608-610: Wzmocnienie Dowodów
**Data:** 2025-12-05
**Cel:** Wzmocnić najsłabsze hipotezy (H7, H1, H9)

---

## Priorytetyzacja Hipotez

### Status obecny:
- ✅ Potwierdzone: 8/12 (H3,H4,H5,H6,H8,H9,H10,H11)
- 🟢 Częściowo: 2/12 (H1,H12)
- 🔴 Niepotwierdzone: 2/12 (**H7**, H2)

### Najsłabsze wymagające wzmocnienia:
1. **H7 (Stałe=Geometria):** 🔴 Niepotwierdzone - QW-601 pokazało tylko ω
2. **H1 (Przestrzeń=Korelacja):** 🟢 Częściowo - d=2.6 nie 3.0
3. **H9 (Grawitacja=Hebbian):** ✅ Ale tylko jakościowe, brak kwantyfikacji

---

## QW-608: Spectral Origin of Constants (H7)

### Hipoteza
α_geo, ω, φ, β pochodzą z **spektralnych właściwości** sieci (eigenvalues).

### Teoria
QW-601 testowało geometryczne kombinacje (π,e,φ,√). Ale może stałe pochodzą z **dynamicznych** właściwości:
- Eigenvalue ratios
- Resonance frequencies
- Spectral gaps

### Test
1. Build coupling matrix S (12×12) z K(d)
2. Compute eigenvalues λᵢ
3. Test ratios:
   - α_geo ≈ λ_max/λ_min?
   - ω ≈ (λᵢ₊₁ - λᵢ)/Δλ_avg?
   - φ ≈ arg(λ_complex)?

**Przewidywanie:** Jeśli stałe = eigenvalue ratios → H7 potwierdzone

**Znaczenie:** Wyprowadzenie α_geo z pierwszych zasad (spektrum sieci)

---

## QW-609: Correlation Dimension Scaling (H1)

### Problem
H1 jest częściowo potwierdzone:
- ✅ Przestrzeń emerguje z korelacji (QW-384, QW-171)
- ❌ Wymiar d=2.6, nie 3.0

### Pytanie
Czy d=3.0 emerguje na **różnych skalach**? (Fractal → smooth transition?)

### Test
Measure correlation dimension d_eff(ρ) w funkcji density/distance:

```python
for rho in [0.01, 0.1, 1.0, 10.0]:  # Różne gęstości
    C(r) = <ψ(x) ψ(x+r)>
    d_eff = -d(log C)/d(log r)
```

**Przewidywanie:**
- Małe ρ (Planck): d≈2.6 (fractal)
- Duże ρ (macro): d→3.0 (smooth limit)

**Znaczenie:** Pokazałoby że d=3 to **emergent effective** na makroskali

---

## QW-610: Multi-Body Hebbian Gravity (H9)

### Problem
H9 potwierdzone tylko **jakościowo** (QW-440: 9.7% kontrakcja).

Brak testów:
- Multi-body (3+ masy)
- Inverse square law (F ∝ 1/r²?)
- Quantitative comparison

### Test
Simulate 3 masses in Hebbian network:

```
Mass A (m=1) at x=0
Mass B (m=1) at x=10  
Mass C (m=0.5) at x=5 (middle)
```

**Measure:**
1. Force on C from A: F_A
2. Force on C from B: F_B
3. Net force: F_net = F_A + F_B

**Test superposition:** F_net = F_A + F_B? (linear!)
**Test inverse square:** F ∝ 1/r²?

**Przewidywanie:** Jeśli Hebbian daje superpozycję + 1/r² → H9 kwantytatywnie potwierdzone

---

## Priorytet Execution

**Najwyższy: QW-608** (H7)
- Najprostsza hipoteza do testowania
- Jasna predykcja (eigenvalue ratios)
- Wypełni największą lukę (H7 całkowicie niepotwierdzone)

**Średni: QW-609** (H1)
- Wzmocni częściowo potwierdzoną H1
- Multi-scale test (elegancki)

**Niski: QW-610** (H9)
- H9 już potwierdzone, to tylko kwantyfikacja
- Mniejszy priorytet niż niepotwierdzone H7

---

**Rekomendacja:** Start z QW-608, potem QW-609.

Gotowy do implementacji?
