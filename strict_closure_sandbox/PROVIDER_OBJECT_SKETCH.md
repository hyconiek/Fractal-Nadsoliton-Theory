## Provider object (Luka 1, N327) – szkic ścisłej konstrukcji

Ten plik jest czysto sandboxowy. Można go bezpiecznie usunąć razem z katalogiem `strict_closure_sandbox`.

### 1. Dane wejściowe z obecnego repo (bez zmiany statusów)

- **Ontologia (AX9/F1)**: nadsoliton jest pojedynczym obiektem informacyjnym.
- **Trasa nad12–sigma (N328–N342)**: dostarcza nośnika `sigma_int` i struktury „residual route” od źródła do obiektów wspierających bridge/export.
- **Luka 1 (N327)**: brakujący obiekt ma być:
  - source-side,
  - observer-free,
  - pair-indexed,
  - noncyclic.

W tym szkicu nie zmieniamy żadnych statusów (`candidate`, `actual`, itd.); tylko proponujemy techniczną formę obiektu zgodną z tym opisem.

### 2. Przestrzeń stanów i indeksów – dwa konkretne modele

Najpierw podajemy ogólny szkielet, a potem dwa konkretne modele referencyjne (do ewentualnego rozwinięcia poza sandboxem).

1. `X` – przestrzeń stanów nadsolitonu:
   - np. przestrzeń mierzalna albo przestrzeń Hilberta, zgodna z dotychczasowymi założeniami o „wewnętrznych stanach”,
   - tutaj pozostawiamy typ przestrzeni jawnie jako parametr (do późniejszej decyzji).
2. `Σ_int` – zbiór (lub przestrzeń) wewnętrznych indeksów `sigma_int`:
   - obiekt przeniesiony z nad12–sigma route,
   - zakładamy jedynie, że istnieje dobrze określona struktura indeksowania (np. skończony lub przeliczalny zbiór, bądź przestrzeń topologiczna).
3. `P_pairs` – zbiór indeksów parowych:
   - np. uporządkowane pary typów (`light`, `matter`), albo bardziej ogólny indeks pary kanałów (abstrakcyjnie: `p ∈ P_pairs`).

Na potrzeby szkicu zakładamy tylko:

- istnienie kartezjańskiego produktu `Σ_int × P_pairs`,
- to, że `X` jest przestrzenią, na której można definiować przekształcenia (operatory, endomorfizmy, funkcje mierzalne itd.).

#### 2.1 Model A (operator liniowy na przestrzeni Hilberta)

- `X` – separowalna przestrzeń Hilberta (np. `L^2` na odpowiedniej przestrzeni konfiguracji nadsolitonu),
- `Σ_int` – przeliczalny zbiór indeksów odpowiadający węzłom trasy nad12–sigma,
- `P_pairs` – skończony zbiór parowych „typów kanałów” (np. `({light}, {matter})`, `({void+light}, {light+matter})` itp.).

W tym modelu każdy `P_{σ,p}` jest operatorem liniowym i ograniczonym na `X`, a dodatkowo możemy wymagać:

- samosprzężoności lub normalności (żeby mieć dobre spektrum),
- kontraktywności (dla kontroli niecykliczności).

#### 2.2 Model B (endomorfizmy przestrzeni mierzalnej)

- `X` – przestrzeń mierzalna z ustaloną miarą (np. `σ`-skończoną),
- `Σ_int` – zbiór indeksów jak wyżej,
- `P_pairs` – jak wyżej,
- `P_{σ,p}` – mierzalne przekształcenia `X → X`, zachowujące (lub kontrakcyjne względem) pewną miarę lub funkcjonał entropijny.

Model B jest bardziej elastyczny, jeśli później chcemy bezpośrednio łączyć providera z argumentami entropijnymi (Shannon-route).

### 3. Definicja: provider jako rodzina operatorów + sanity check

**Definicja (szkicowa).**  
Provider object jest dany przez rodzinę operatorów

```text
  { P_{σ,p} : X → X }_{(σ,p) ∈ Σ_int × P_pairs}
```

spełniających następujące warunki:

1. **Source-side**  
   - każdy `P_{σ,p}` jest zdefiniowany niezależnie od obserwatora; dziedzina i przeciwdziedzina to `X` (wewnętrzne stany nadsolitonu),
   - żadne parametry obserwacyjne (zewnętrzne dane, wybór układu odniesienia) nie wchodzą do definicji `P_{σ,p}`.

2. **Observer-free**  
   - dla dowolnej transformacji obserwatora `O` działającej jako automorfizm na `X` (czyli `O : X → X` odzwierciedlający zmianę „opisu” stanu, nie samego stanu),
   - wymagamy
     ```text
     O ∘ P_{σ,p} = P_{σ,p} ∘ O
     ```
     lub, słabiej, że klasy równoważności stanów pod działaniem obserwatora są zachowywane przez `P_{σ,p}`.

3. **Pair-indexed**  
   - indeks `p ∈ P_pairs` jest nieusuwalną częścią obiektu:
     - dla różnych `p1 ≠ p2` rodziny `{P_{σ,p1}}_σ` i `{P_{σ,p2}}_σ` są rozróżnialne (np. różnią się spektrum, strukturą orbit albo typem generowanych obiektów),
   - w praktyce:
     - można to zrealizować jako różne klasy operatorów (np. „światło-dominujące” vs „materia-dominujące”), ale ten szczegół pozostaje na późniejszy etap.

4. **Noncyclic**  
   - dla każdej pary `(σ,p)` oraz prawie każdego (w sensie miary / topologii) stanu `x ∈ X`:
     ```text
     P_{σ,p}^n(x) ≠ x  dla każdego n ≥ 1
     ```
     lub, słabiej, cykle są mierzalnie/rzadkie i nie dominują dynamiki,
   - w sensie abstrakcyjnym: provider generuje trajektorie, które „idą naprzód” w strukturze stanów, zamiast wracać do punktu wyjścia w skończonej liczbie kroków.

Warunki 1–4 odpowiadają odpowiednio: source-side, observer-free, pair-indexed, noncyclic.

#### 3.1 Sprawdzenie spójności warunków w Modelu A

W Modelu A (Hilbertowskim) możemy wybrać:

- klasę dopuszczalnych transformacji obserwatora `O` jako grupę unitarnych operatorów na `X`,
- providera jako rodzinę operatorów, które:
  - albo komutują z tą grupą (silne observer-free),
  - albo należą do jej komutantu.

**Obserwacja:**  
Jeśli `P_{σ,p}` leży w komutancie grupy unitarnych obserwatorów, to warunek observer-free jest spełniony w sensie algebraicznym:

```text
  O P_{σ,p} = P_{σ,p} O  dla każdego obserwatora O
```

Noncykliczność można uzyskać np. przez:

- wybór operatorów o spektum zawartym w tarczy jednostkowej z wyjątkiem punktów pierwiastków z jedności,
- dodanie lekkiej kontrakcji (np. norma operatora < 1), co eliminuje cykle skończone poza zdegenerowanymi przypadkami.

#### 3.2 Sprawdzenie spójności w Modelu B

W Modelu B (przestrzeń mierzalna) observer-free można rozumieć jako:

- grupa transformacji obserwatora `O` działa jako automorfizmy mierzalne zachowujące klasę miary,
- `P_{σ,p}` zachowuje klasy równoważności prawie wszędzie względem tej miary (czyli trajektorie różniące się tylko „przez obserwatora” pozostają równoważne).

Noncykliczność można w tym modelu oprzeć na:

- wymaganiu, że pewien funkcjonał (np. entropia lub „głębokość fraktalna”) jest ściśle rosnący wzdłuż orbit `P_{σ,p}`,
- co automatycznie wyklucza cykle skończone poza przypadkami miary zero.

### 4. Jak to łączy się z nad12–sigma route

Zamiast modyfikować istniejące pakiety, zakładamy następujące relacje (do późniejszego udowodnienia w głównym repo):

1. **Zgodność z `sigma_int`**  
   - dla każdego `σ ∈ Σ_int` trasa nad12–sigma określa pewną klasę stanów lub konfiguracji w `X`,
   - operator `P_{σ,p}` powinien działać „wzdłuż” tej klasy, tj.:
     - zachowywać przynależność do odpowiedniego ramienia trasy,
     - albo przenosić stany między węzłami trasy zgodnie z już istniejącą strukturą.

2. **Zasilanie bridge/export**  
   - trajektorie generowane przez iteracje `P_{σ,p}`:
     ```text
     x, P_{σ,p}(x), P_{σ,p}^2(x), ...
     ```
     powinny dostarczać kandydatów na:
     - obiekty wspierające residual bridge (Luka 2),
     - populacje parowe potrzebne w dalszych krokach (theta/pair-population candidates).

Formalny krok w głównym repo polegałby na:

- zdefiniowaniu konkretnej klasy operatorów `P_{σ,p}` (np. jako operatorów liniowych/bounded na `X`, albo jako transformacji mierzalnych),
- udowodnieniu, że:
  - są one zgodne z istniejącą nad12–sigma infrastrukturą,
  - generują strukturę wymaganą przez N302 itd.

W Modelu A/B można to precyzyjniej sformułować np. jako:

- istnienie odwzorowania:
  ```text
    Φ : Σ_int × P_pairs → Klasa_tras_nad12−sigma
  ```
  takie, że orbity `P_{σ,p}` są zanurzone (immersed) w odpowiadających im trasach z obrazu `Φ`.

### 5. Minimalny pakiet twierdzeń docelowych

Docelowo (po przeniesieniu poza sandbox) ścieżka mogłaby wyglądać tak:

1. **Definicja formalna**  
   - pełna definicja `X`, `Σ_int`, `P_pairs`, klasy operatorów `P_{σ,p}`.
2. **Twierdzenie 1 (source-side / observer-free)**  
   - dowód, że dla dowolnych dopuszczalnych transformacji obserwatora `O` zachodzi komutacja (lub zachowanie klas równoważności), czyli brak zależności od obserwatora.
3. **Twierdzenie 2 (pair-indexed)**  
   - dowód, że różne `p` dają nieredukowalne (w sensie nieizomorficzne) klasy operatorów lub trajektorii.
4. **Twierdzenie 3 (noncyclic)**  
   - dowód, że prawie żadna trajektoria nie jest skończenie cykliczna (z dokładnie określonym sensem „prawie”).
5. **Twierdzenie 4 (bridge-feeding)**  
   - dowód, że trajektorie `P_{σ,p}` dostarczają obiektów, które spełniają warunki wejściowe do konstrukcji residual bridge/export (Luka 2).

Ten plik nie rozwiązuje tych twierdzeń – tylko ustawia ramę, w której można je formułować i dowodzić bez naruszania istniejących guardraili (brak cichego mostu legacy↔strict, brak ukrytego selektora).

---

### 6. Prototyp konkretnego providera w Modelu A

Poniżej podajemy jeden **prosty, ale niebanalny** przykład providera w Modelu A, który spełnia warunki 1–4 w sensie matematycznym. To jest tylko prototyp – nie jest jeszcze zszyty z pełną semantyką nad12–sigma, ale pokazuje, że klasa konstrukcji nie jest pusta.

#### 6.1 Ustawienie

- `X = ℓ^2(ℤ)` – przestrzeń Hilberta wszystkich ciągów kwadratowosumywalnych `ψ = (ψ_k)_{k∈ℤ}`,
- obserwatorzy: grupa unitarnych operatorów przesunięcia fazy
  ```text
    (O_θ ψ)_k = e^{i θ} ψ_k,   θ ∈ [0, 2π)
  ```
  oraz (opcjonalnie) przesunięcia indeksu
  ```text
    (S ψ)_k = ψ_{k-1}.
  ```
- `Σ_int` – dowolny przeliczalny zbiór, który w tym prototypie interpretujemy jako etykiety „ramion” (ale nie wchodzimy w szczegóły),
- `P_pairs` – skończony zbiór, np.
  ```text
    P_pairs = {p_LM, p_ML}
  ```
  dla „light→matter” i „matter→light”.

#### 6.2 Definicja operatorów `P_{σ,p}`

Wybieramy jedno ramię `σ ∈ Σ_int` i definiujemy (dla przejrzystości) operator niezależny od `σ`, ale różny dla różnych `p`:

1. **Dla `p = p_LM`** (kanał „light→matter”):
   ```text
   (P_{σ,p_LM} ψ)_k = a ψ_{k-1},   gdzie 0 < |a| < 1.
   ```
   Jest to „przesunięcie z tłumieniem” w prawo.

2. **Dla `p = p_ML`** (kanał „matter→light”):
   ```text
   (P_{σ,p_ML} ψ)_k = b ψ_{k+1},   gdzie 0 < |b| < 1.
   ```
   Jest to „przesunięcie z tłumieniem” w lewo.

Można też dopuścić zależność `a = a(σ)`, `b = b(σ)`, ale nie jest to konieczne do spełnienia warunków 1–4.

#### 6.3 Weryfikacja warunków 1–4

1. **Source-side**  
   - `P_{σ,p}` działa wyłącznie na `X = ℓ^2(ℤ)`, nie używa żadnych zewnętrznych parametrów obserwacyjnych.

2. **Observer-free**  
   - dla obserwatorów fazowych `O_θ` mamy:
     ```text
       (O_θ P_{σ,p_LM} ψ)_k = e^{i θ} a ψ_{k-1},
       (P_{σ,p_LM} O_θ ψ)_k = a e^{i θ} ψ_{k-1},
     ```
     czyli `O_θ P_{σ,p_LM} = P_{σ,p_LM} O_θ`, analogicznie dla `p_ML`;  
   - jeśli chcemy uwzględnić również przesunięcia indeksu `S`, możemy zawęzić klasę obserwatorów do samych `O_θ`, albo lekko zmodyfikować definicję `P_{σ,p}` tak, by należały do komutantu wybranej grupy obserwatorów.

3. **Pair-indexed**  
   - operatory dla `p_LM` i `p_ML` są różne (działają w przeciwnych kierunkach, z ewentualnie różnymi współczynnikami `a`, `b`), więc rodziny `{P_{σ,p_LM}}_σ` i `{P_{σ,p_ML}}_σ` są nieredukowalne do jednej klasy poprzez automorfizmy zachowujące orientację indeksów;  
   - w szczególności, spektra i własności ergodyczne mogą się różnić dla różnych `p`.

4. **Noncyclic**  
   - ponieważ `0 < |a|, |b| < 1`, operator `P_{σ,p}` jest ściśle kontrakcyjny w normie Hilberta:
     ```text
       ∥P_{σ,p} ψ∥ = |a| ∥ψ∥  lub  |b| ∥ψ∥  < ∥ψ∥  dla ψ ≠ 0.
     ```
   - jeśli istniałby cykl skończony `P_{σ,p}^n ψ = ψ` z `ψ ≠ 0`, to mielibyśmy:
     ```text
       ∥ψ∥ = ∥P_{σ,p}^n ψ∥ = |a|^n ∥ψ∥  (lub |b|^n ∥ψ∥),
     ```
     co implikuje `|a|^n = 1` (lub `|b|^n = 1`), a to jest sprzeczne z `0 < |a|,|b| < 1`;  
   - zatem **żaden** niezerowy stan nie może leżeć na cyklu skończonym – noncykliczność zachodzi globalnie (nie tylko „prawie wszędzie”).

Ten prototyp pokazuje, że:

- istnieją rodziny operatorów `P_{σ,p}` spełniające 1–4 dokładnie w modelu operatorowym,
- observer-free i noncykliczność można uzyskać zwykłymi technikami z analizy operatorowej (komutant grupy unitarnej + kontrakcja),
- pair-indexing może być nośny już na poziomie prostych różnic w kierunku działania i parametrach tłumienia.

Kolejny krok (poza sandboxem) to:

- dopasowanie indeksu `k` i tłumienia `a(σ), b(σ)` do realnej struktury węzłów nad12–sigma,
- zdefiniowanie jawnego odwzorowania `Φ(σ, p)` na konkretne trasy z N328–N342, tak aby trajektorie `P_{σ,p}` miały sens fizyczny/ontologiczny zgodny z resztą FAR.

