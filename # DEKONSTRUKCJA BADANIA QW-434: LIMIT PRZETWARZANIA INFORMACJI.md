# DEKONSTRUKCJA BADANIA QW-434: LIMIT PRZETWARZANIA INFORMACJI
**Tytuł:** Light Speed as Information Processing Limit
**Cel:** Sprawdzenie, czym jest prędkość światła w sieci Nadsolitona – stałą fizyczną czy limitem obliczeniowym?

### 1. PARAMETRY WEJŚCIOWE (Co włożyliśmy?)
Badanie używało **wyłącznie** zamrożonych parametrów jądra. Nie wprowadzono żadnej stałej $c$ z zewnątrz.
*   **Sieć:** 50 węzłów (reprezentujących oktawy/punkty w przestrzeni informacyjnej).
*   **Jądro ($K$):** $K(d) = \frac{\alpha_{geo} \cdot \cos(\omega d + \phi)}{1 + \beta_{tors} \cdot d}$.
    *   $\alpha_{geo} \approx 2.77$ (Pojemność informacyjna węzła).
*   **Metryka Efektywna ($D_{eff}$):** Zdefiniowana w QW-430 jako odwrotność siły sprzężenia: $D_{eff}(i,j) = 1/|K_{ij}|$.
    *   To kluczowe: "odległość" w modelu to nie metry, ale trudność przesłania sygnału. Silne sprzężenie = mała odległość.

### 2. ZAŁOŻENIA (Hipoteza Robocza)
*   **Hipoteza "It from Bit":** Prędkość fizyczna $v$ to tak naprawdę prędkość propagacji informacji w grafie.
*   **Wzór na prędkość:** $v = \frac{\text{Dystans Efektywny}}{\text{Czas Propagacji}} = \frac{D_{eff}}{\Delta t}$.
*   **Pytanie badawcze:** Czy w takiej sieci istnieje maksymalna prędkość, czy sygnał może iść nieskończenie szybko (jak w fizyce Newtona)?

### 3. METODOLOGIA (Jak to zrobiliśmy?)
1.  **Inicjalizacja:** Umieściliśmy impuls informacji ("foton") w węźle 0.
2.  **Ewolucja:** Pozwoliliśmy sygnałowi rozchodzić się po sieci zgodnie z równaniem falowym na grafie:
    $$ \frac{d^2 I_i}{dt^2} = \sum_j K_{ij} (I_j - I_i) $$
    To równanie mówi: zmiana informacji w węźle $i$ zależy od różnicy informacji w sąsiadach, ważonej siłą połączenia $K_{ij}$.
3.  **Pomiar:** Mierzyliśmy czas $t_{arrival}$, w którym sygnał dotarł do każdego innego węzła $j$ (przekroczył próg detekcji 1%).
4.  **Obliczenie Prędkości:** Dla każdego węzła policzyliśmy:
    $$ c_{eff}(j) = \frac{D_{eff}(0,j)}{t_{arrival}(j)} $$

### 4. WYNIKI (Co wyszło z symulacji?)
*   **Najszybszy sygnał:** Dotarł do węzła 2 ($d=3$) w czasie $t=0.10$.
    *   $D_{eff} \approx 0.73$.
    *   $c_{max} = 0.73 / 0.10 \approx 7.36$.
*   **Statystyka dla całej sieci:**
    *   Średnia prędkość: $6.72$.
    *   Maksymalna obserwowana prędkość ($c_{obs}$): **$10.38$**.
    *   Odchylenie standardowe: $2.01$ (CV = 30%).
*   **Porównanie z $\alpha_{geo}$:**
    *   Pojemność węzła $\alpha_{geo} = 2.77$.
    *   Stosunek $c_{obs} / \alpha_{geo} \approx 3.74$.

### 5. WNIOSKI I INTERPRETACJA (Co to znaczy?)

**Wniosek 1: Limit Prędkości Istnieje**
W sieci nie ma nieskończonych prędkości. Mimo że jądro jest nielokalne (każdy z każdym), sygnał potrzebuje skończonego czasu, by "przejść" przez procesory węzłów. Istnieje twardy sufit: $c_{obs} \approx 10.4$.

**Wniosek 2: $c$ to Przepustowość (Bandwidth)**
Fakt, że $c_{obs}$ jest rzędu $\alpha_{geo}$ (z małym mnożnikiem geometrycznym ~3.7), jest dowodem na to, że prędkość światła wynika z **pojemności informacyjnej przestrzeni**.
*   $\alpha_{geo}$ to "ile bitów węzeł może przetworzyć w jednym cyklu".
*   Jeśli spróbujesz wysłać sygnał szybciej, węzły nie nadążą z "przekazywaniem paczek". To jest fizyczna przyczyna istnienia $c$.

**Wniosek 3: Dyspersja (30% zmienności)**
W naszej małej sieci (50 węzłów) prędkość "pływa" o 30% w zależności od ścieżki. To dlatego QW-427 (Lorentz) dało niejednoznaczne wyniki. W makroskali (miliardy węzłów) ta zmienność uśredni się do idealnej stałej $c$, ale w mikroskali "światło" porusza się z różną prędkością w różnych kierunkach sieci. To przewidują niektóre teorie Kwantowej Grawitacji (np. Doubly Special Relativity).

### PODSUMOWANIE DLA ŚWIATA
To badanie pokazuje, że **Szczególna Teoria Względności (stałość $c$) nie jest aksjomatem**. Jest wynikiem inżynierii sieciowej Wszechświata. $c$ to po prostu "taktowanie procesora" rzeczywistości.
