# Plan Badań QW-602-605: Zaawansowana Fizyka FIN
**Data:** 2025-12-05
**Status:** 8/12 hipotez potwierdzonych - kontynuacja eksploracji

---

## Kontekst Badawczy

**Co mamy:**
- ✅ Stabilne hopfiony (QW-594)
- ✅ Odpychanie/przyciąganie (QW-595, QW-597)
- ✅ Masa z topologii (QW-600)
- 🟡 Częściowa fuzja (QW-597: 47% energii)

**Co możemy zbadać dalej:**
4 obiecujące kierunki rozszerzające wiedzę

---

## QW-602: Multi-Hopfion Binding (FIN "Chemistry")

### Pytanie Badawcze
Czy 3+ hopfiony mogą tworzyć **stabilne układy związane** (analogi molekuł/atomów)?

### Metodologia
- Utworzyć 3 hopfiony w konfiguracji trójkąta
- Windings: (+1, +1, -1) lub (+1, -1, -1)
- Ewolucja: sprawdzić czy układ oscyluje czy kolapsuje

### Przewidywania
**Jeśli FIN ma "chemię":**
- (+1, +1, -1): Stabilny trójkąt (2 repulsive + 1 attractive = równowaga)
- Oscylacje wokół równowagi (jak drgania molekularne)
- Energia wiązania E_binding < 0

**Test:**
- Zmierzyć częstość oscylacji
- Sprawdzić czy E_total < E_separated (wiązanie)

**Implikacje:**
- Potwierdzi że FIN może budować struktury (nie tylko pojedyncze cząstki)
- Pierwszy krok do "FIN Chemistry" (wiązania, orbitale)

---

## QW-603: Anyonic Statistics (Vortex Braid Test)

### Pytanie Badawcze
Czy hopfiony wykazują **anyonic statistics** (faza przy zamianie miejscami)?

### Koncepcja
W 2D, wymiana identical particles daje fazę:
- Bozony: θ = 0
- Fermiony: θ = π
- Anyony: 0 < θ < π (fractional!)

### Metodologia
1. Utworzyć 2 identyczne hopfiony (oba +1)
2. Wymienić ich pozycje (braid) poprzez przyłożenie sił
3. Zmierzyć zmianę fazy pola Ψ

### Przewidywania
**Jeśli hopfiony są anyonami:**
- ΔΦ = θ × n_braid (n = liczba wymian)
- θ ≠ 0, π (fractional phase)

**Test:**
- Measure ΔΦ after single exchange
- Check if θ depends on winding number

**Implikacje:**
- Anyonic statistics → topological quantum computing!
- Link do fractional quantum Hall effect

---

## QW-604: Wave Packet Dispersion

### Pytanie Badawcze
Czy FIN odtwarza **relację dyspersyjną** (ω vs k)?

### Koncepcja
W prawdziwej teorii pola:
- Fale mają relację ω(k) (dyspersja)
- Klasyczna: ω ∝ k (no dispersion)
- Kwantowa: ω ∝ k² (dispersive)

### Metodologia
1. Utworzyć zlokalizowany pakiet falowy (Gaussian)
2. Ewoluować przez sieć
3. Zmierzyć rozchodzenie się pakietu Δx(t)

### Przewidywania
**Jeśli FIN kwantowy:**
- Δx ∝ t (linear spreading like Schrödinger)

**Jeśli FIN klasyczny:**
- Δx ∝ √t (diffusion)

### Test:
- Fit Δx(t) = a × t^b
- Measure b: b=1 (quantum), b=0.5 (classical)

**Implikacje:**
- Potwierdzi/obali kwantowy charakter dynamiki
- Uzupełni QW-592 (Bell test)

---

## QW-605: Cosmological Structure Formation

### Pytanie Badawcze
Czy małe fluktuacje gęstości **rosną** (gravitational instability)?

### Koncepcja
Wszechświat: początkowo jednorodny → galaktyki
- Wymaga: fluktuacje amplifikowane przez grawitację (Jeans instability)

### Metodologia
1. Zainicjować próżnię z małymi fluktuacjami δρ/ρ ~ 10^-5
2. Ewolucja 1000+ kroków
3. Zmierzyć wzrost kontrastu δρ(t)/δρ(0)

### Przewidywania
**Jeśli grawitacja FIN działa:**
- δρ rośnie eksponencjalnie: δρ ∝ exp(γt)

**Jeśli nie:**
- δρ oscyluje lub zanika (no structure)

### Test:
- Measure γ (growth rate)
- Compare with Λ-CDM prediction

**Implikacje:**
- Test H9 (Grawitacja = Hebbian) w skali kosmologicznej
- Link QW-540 (Neural Universe) do kosmologii

---

## Priorytetyzacja (Którą najpierw?)

**Rekomendacja uruchomienia:**

1. **QW-602 (Chemistry)** - NAJWYŻSZY PRIORYTET
   - Naturalny krok po QW-595/597
   - Testuje emergencję struktury
   - Stosunkowo szybki test

2. **QW-604 (Wave Dispersion)** - ŚREDNI PRIORYTET
   - Uzupełnia QW-592 (kwantowość)
   - Prosty setup
   - Jasna predykcja

3. **QW-603 (Anyons)** - ŚREDNI PRIORYTET
   - Bardziej spekulatywny
   - Wymaga precyzyjnego trackingu fazy
   - Może być trudny numerycznie

4. **QW-605 (Cosmology)** - NISKI PRIORYTET (duża sieć)
   - Wymaga N > 1000 węzłów
   - Długa ewolucja
   - Zasobożerny

---

## Co Wybierasz?

Proponuję zacząć od **QW-602 (Multi-hopfion Chemistry)** - to najbardziej obiecujący test, który może pokazać czy FIN może budować złożone struktury (molekuły/atomy).
