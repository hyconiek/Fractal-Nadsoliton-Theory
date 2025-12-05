# Plan Badawczy: Droga Naprzód - Kwantyzacja i Rekurencja
**Cel:** Przekształcenie teorii z modelu deterministycznego (ODE) na w pełni kwantowy (Path Integrals) i fraktalny (Recursive Nesting).

## Filary Nowego Podejścia

### 1. Kwantyzacja: Od ODE do Path Integrals
*   **Problem:** Dotychczasowe symulacje ($dA/dt$) badały jedną, "klasyczną" trajektorię.
*   **Rozwiązanie:** Sumowanie po wszystkich możliwych ścieżkach w grafie. Amplituda przejścia $K(a, b) = \sum_{paths} e^{i S[x(t)]}$.
*   **Hipoteza:** Klasyczna geodezyjna (orbita) wyłoni się jako linia najmniejszego działania (konstruktywna interferencja).

### 2. Pełne Zagnieżdżenie (Recursive Nesting)
*   **Problem:** Symulowaliśmy warstwy osobno ("Warstwa 0", "Warstwa 1").
*   **Rozwiązanie:** Symulacja, w której węzeł sieci *jest* całą pod-siecią.
*   **Implementacja:** Obiekt `FractalNode` zawierający w sobie instancję `Simulation`.

### 3. Statystyczna Mechanika Grawitacji ($1/r^2$)
*   **Problem:** Grawitacja była modelowana jako przepływ płynu (Navier-Stokes).
*   **Rozwiązanie:** Wyprowadzenie $1/r^2$ z czystej statystyki błądzenia losowego na grafie fraktalnym.
*   **Hipoteza:** Zagęszczenie ścieżek w pobliżu masy powoduje efektywne przyciąganie (siła entropowa).

---

## Planowane Eksperymenty (Jupyter Notebook / Colab)

### **QW-580: Path Integral Propagator**
*   Obliczenie propagatora $K(x, t)$ na sieci dyskretnej.
*   Wizualizacja interferencji fal prawdopodobieństwa.

### **QW-581: Emergent Gravity form Paths**
*   Pomiar prawdopodobieństwa znalezienia cząstki próbnej w pobliżu masy, gdy sumujemy miliony ścieżek losowych.
*   Sprawdzenie czy $P(r) \sim 1/r$ (potencjał).

### **QW-582: Recursive Simulation (Toy Model)**
*   Symulacja 2-poziomowa: Sieć "Wszechświata", gdzie każdy węzeł to "Galaktyka" (też sieć).
*   Demonstracja przepływu czasu: Czy czas w "Galaktyce" płynie szybciej niż we "Wszechświecie"?

## Instrukcja dla Użytkownika (Colab)
Ze względu na złożoność obliczeniową (sumowanie ścieżek to $O(N^k)$ lub Monte Carlo), kod zostanie przygotowany jako **Jupyter Notebook** do uruchomienia w Google Colab.
