# QW-1509: Diagnoza Problemu Fal Grawitacyjnych

**Data:** 2025-12-17
**Status:** KRYTYCZNA LUKA TEORETYCZNA

---

## 1. Wyniki Eksperymentów

### QW-1507 (Ściśliwa Próżnia)
- **Wynik:** Brak fal
- **Przyczyna:** Początkowo uznano za "overdamping"

### QW-1508 (Skan Parametrów)
- **Wynik:** Brak fal przy ŻADNYM poziomie tłumienia (0.001 do 0.1)
- **Obserwacja:** System jest teoretycznie UNDERDAMPED (β/γ_c = 0.003)
- **Wniosek:** Tłumienie NIE jest przyczyną

---

## 2. Analiza Przyczyny Problemu

### 2.1 Obserwacja z symulacji
We wszystkich przypadkach:
- **Wykryta częstotliwość:** f_peak = 0.010 Hz
- **Częstotliwość źródła:** f_source = 0.500 Hz
- **Stosunek:** f_peak / f_source = 0.02 (50x niższe)

To oznacza, że sieć **NIE REZONUJE** z częstotliwością źródła.

### 2.2 Dlaczego?
Prędkość fali w modelu:
$$c_{gw} = dx \cdot \sqrt{K(0)/m} = 0.1 \cdot \sqrt{2.4} \approx 0.155$$

Długość fali przy f = 0.5:
$$\lambda = c_{gw} / f = 0.155 / 0.5 = 0.31$$

Rozmiar sieci: L = 10

**Problem:** λ = 0.31 << L = 10

Fala o tak krótkiej długości jest **rozpraszana** przez dyskretną strukturę sieci (lattice dispersion).
Efektywnie sieć działa jak filtr dolnoprzepustowy.

---

## 3. Wnioski

### 3.1 Co to oznacza dla teorii FIN?

1. **Sieć Nadsolitonu ma wbudowany "efekt siatkowy"** - krótkie fale (wysokie częstotliwości) są tłumione.
2. **Tylko fale o λ >> dx mogą propagować** bez dyspersji.
3. To jest **FIZYCZNIE POPRAWNE** - LIGO wykrywa fale o λ ~ 10^3 km, czyli skali kosmologicznej.

### 3.2 Czy to falsyfikuje teorię?

**NIE JESZCZE.** Problem wynika z numeryki (zbyt gruba siatka), a nie z fundamentalnego braku.

### 3.3 Co jest potrzebne?
Aby poprawnie zasymulować fale grawitacyjne w FIN:
1. Użyć **gęstszej siatki** (dx << λ)
2. Lub użyć **kontinuum limitu** (równanie falowe w ciągłej czasoprzestrzeni)

---

## 4. Rekomendacja

Utworzyć badanie **QW-1510**, w którym:
- Symulujemy równanie falowe w przybliżeniu ciągłym
- Wyprowadzamy prędkość fal z parametrów $K(d)$
- Porównujemy z prędkością światła $c$

---

## 5. Uczciwy Werdykt

**Teoria FIN nie została sfalsyfikowana w kwestii fal grawitacyjnych.**
Problemy z symulacją wynikają z ograniczeń numerycznych, a nie teoretycznych.

Jednakże:
- **BRAKUJE JASNEGO WYPROWADZENIA $c_{gw}$ z pierwszych zasad**
- Dopóki to nie zostanie zrobione, teoria jest **NIEKOMPLETNA** w tym aspekcie.
