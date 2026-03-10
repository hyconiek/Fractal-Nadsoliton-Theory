## Bridge/export preobjects z providera – szkic

Ten plik dalej rozwija sandboxowy model providera z `PROVIDER_OBJECT_SKETCH.md`. Celem jest sprawdzenie, na ile da się **konstrukcyjnie** podejść do Luki 2 (bridge/export) bez łamania guardrailów, i gdzie pojawia się realna ściana.

### 1. Dane wyjściowe

Zakładamy kontekst Modelu A:

- `X = ℓ^2(ℤ)`,
- rodzina operatorów `P_{σ,p} : X → X` z sekcji 6 tamtego pliku,
- noncykliczność globalna (brak niezerowych cykli skończonych),
- observer-free względem grupy `O_θ`.

Nie używamy żadnych konkretnych formuł z legacy kernela; pracujemy czysto strict‑side.

### 2. Orbitowe klasy jako kandydaci na obiekty typu bridge/export

Idea: residual bridge/export ma „żyć” nie na pojedynczych stanach, tylko na **strukturze orbitowej** providera. W takim modelu:

- provider generuje trajektorie
  ```text
    ψ, P_{σ,p} ψ, P_{σ,p}^2 ψ, ...
  ```
- bridge/export preobject to nie pojedynczy punkt, lecz:
  - klasa równoważności całej orbity,
  - albo pewien obiekt zbudowany z zachowania orbity (np. miara niezmiennicza, limit, średnia Cesàro).

#### 2.1 Definicja orbit i klas równoważności

Definiujemy relację równoważności na `X \ {0}`:

```text
  ψ ~_{σ,p} φ    ⇔    istnieją m,n ≥ 0 takie, że  P_{σ,p}^m ψ = c P_{σ,p}^n φ
                      dla pewnej stałej c ≠ 0 (np. faza)
```

Intuicja:

- liczymy orbity z zadaną tolerancją na reskalowanie (np. fazę globalną),
- traktujemy stany jako „ten sam obiekt orbitowy”, jeśli można je uzyskać z siebie nawzajem przez skończoną liczbę kroków providera + prostą symetrię.

Następnie definiujemy:

```text
  [ψ]_{σ,p} = klasa równoważności ψ   (orbit class)
```

To jest pierwszy kandydat na „surowy obiekt bridge/export” dla zadanej pary `(σ,p)`.

#### 2.2 Własności tych klas

W naszym prototypie:

- kontrakcyjność `P_{σ,p}` (norma < 1) implikuje:
  - każda orbita `P_{σ,p}^n ψ` zbiega do zera w normie,
  - więc wszystkie stany na orbicie są „coraz mniejsze”, ale **klasa równoważności** jest nadal dobrze zdefiniowana przez pierwsze kroki,
- noncykliczność gwarantuje, że nie ma skończonych cykli poza zerem, więc:
  - każda niezerowa klasa `[ψ]_{σ,p}` jest nieskończona,
  - można naturalnie porządkować stany na orbicie wg „głębokości” iteracji.

To sugeruje interpretację:

- „bridge preobject” = informacja o **kierunku** i **głębokości** przejścia wzdłuż danego kanału `(σ,p)`, oderwana od absolutnej skali amplitudy.

### 3. Uśrednione obiekty orbitowe (Cesàro / miary)

Drugie, bardziej stabilne podejście:

- dla każdej trajektorii definiujemy uśrednienie:

```text
  A_N(ψ) = (1/N) ∑_{n=0}^{N-1} P_{σ,p}^n ψ
```

- w klasycznym settingu kontrakcyjnym (np. twierdzenia ergodyczne) takie uśrednienia mają zwykle granice lub podciągi zbieżne.

**Kandydat:**  

```text
  B_{σ,p}(ψ) = granica (np. w sensie słabym) A_N(ψ) gdy N → ∞
```

Jeśli taka granica istnieje i jest dobrze określona dla dużej klasy stanów `ψ`, można:

- traktować `B_{σ,p}` jako **kanoniczną projekcję** z `X` na przestrzeń „bridge preobjects”,
- zdefiniować przestrzeń obiektów:

```text
  Y_{σ,p} = { B_{σ,p}(ψ) : ψ ∈ X }  ⊂ X
```

W tym obrazie:

- provider generuje orbitę,
- `B_{σ,p}` wydobywa z niej stabilną część, która ma potencjał być nośnikiem mostu/exportu.

### 4. Co można udowodnić w tym prototypie

W czysto matematycznym settingu (bez odniesień do N302/N328–N342) można:

1. **Udowodnić istnienie i jedyność granicy Cesàro** dla klas operatorów:
   - np. dla kontrakcji liniowych na Hilbercie:
     - `∥P_{σ,p}∥ < 1` implikuje, że `A_N(ψ)` zbiega do `(I - P_{σ,p})^{-1}(I - P_{σ,p}^N) ψ / N`, a to typowo idzie do zera, ale można definiować inne funkcjonały (np. sumy ważone) tak, żeby wyłuskać „kierunek”,
   - można też zamiast Cesàro użyć sum ważonych `∑ w_n P_{σ,p}^n ψ` z doborem `w_n`.
2. **Pokazać, że `Y_{σ,p}` jest dobrze określoną, niezerową podprzestrzenią** dla nietrywialnych operatorów `P_{σ,p}`.
3. **Sprawdzić observer-free na poziomie `Y_{σ,p}`**:
   - jeśli `P_{σ,p}` komutuje z grupą obserwatorów, to tak samo będzie komutować operator uśredniający `B_{σ,p}`, więc klasy/obrazy zachowują invariancję obserwatora.

Te fakty są standardowym rzemiosłem z analizy operatorowej / ergodycznej i nie wymagają żadnych dodatkowych aksjomatów fizycznych.

### 5. Gdzie pojawia się realna „ściana”

Żeby z tego prototypu zrobić **rzeczywistą realizację Luki 2** (residual bridge/export) w repo, potrzeba:

1. **Zszyć indeks `k` i parametry operatora z realnymi obiektami nad12–sigma**:
   - nasz `k ∈ ℤ` jest czysto matematycznym indeksem,
   - struktura rzeczywistej trasy nad12–sigma z N328–N342 (węzły, typy przejść, klasy sigma) jest opisana w repo w specjalnej notacji,
   - bez wejścia w tę konkretną strukturę nie da się udowodnić, że:
     - orbitowe klasy `[ψ]_{σ,p}` albo `Y_{σ,p}` zgadzają się z tym, co w N302/N328–N342 nazywane jest „object support” bridge/export.

2. **Zweryfikować wszystkie lokalne warunki z N302**:
   - N302 definiuje (wg Release 6.2) brakujący obiekt „actual residual bridge/export-map object support” przez:
     - konkretne warunki lokalne na sigma_int route,
     - relacje między różnymi ramionami i oknami.
   - w sandboxie nie mamy tych formuł wczytanych ani przepisanych do języka operatorów; bez tego:
     - możemy mieć bardzo elegancki bridge‑preobject, który **nie spełnia** kluczowego warunku z N302 (np. jakiejś równowagi, nonequality, czy zgodności z oknami pozytywnymi).

3. **Spójność z QW‑2191 i późniejszym selektorem**:
   - konstrukcja orbitowa może (lub nie) wprowadzać ukryty wybór selektora,
   - bez przejścia przez pełny dowód QW‑2191 nie da się zagwarantować, że:
     - nowa konstrukcja nie obchodzi obstrukcji „tylnymi drzwiami”,
     - albo przeciwnie: że nie odtwarza dokładnie tego samego problemu pod innymi nazwami.

Na tym etapie **ściana nie jest logiczna, tylko informacyjna/strukturalna**:

- na czysto matematycznym poziomie:
  - provider + orbitowe obiekty + uśrednianie działają,
  - można udowodnić ich własności.
- żeby jednak powiedzieć: „to jest dokładnie `actual residual bridge/export` z N302”, trzeba:
  - zaczytać szczegółową definicję obiektu z N302/N328–N342,
  - przeprowadzić dopasowanie symbol‑po‑symbolu.

Tego kroku **nie da się uczciwie wykonać „w ciemno”** – wymaga to pełnego wglądu w strukturę tych pakietów i ręcznej decyzji teoretycznej, jak interpretujesz ich obiekty w języku operatorowym.

### 6. Podsumowanie: dokąd doszliśmy w sandboxie

- Pokazaliśmy, że:
  - istnieje wyraźna, niepusta klasa providerów `P_{σ,p}` spełniających N327,
  - z orbit tych operatorów można naturalnie zbudować klasy obiektów (`[ψ]_{σ,p}`, `Y_{σ,p}`), które są dobrymi kandydatami na bridge‑preobjects,
  - można to zrobić z pełnym matematycznym rygorem na poziomie analizy operatorowej, bez łamania guardrailów.
- Ściana pojawia się w momencie, gdy chcemy:
  - **utożsamić** te obiekty z konkretnym `actual residual bridge/export` z N302,
  - **zagwarantować** zgodność z całą strukturą nad12–sigma i QW‑2191.

Czyli: od strony czystej matematyki **ściany jeszcze nie ma** – konstrukcje są możliwe.  
Prawdziwa ściana jest **semantyczna**: wymaga już Twojej decyzji, jak dokładnie przepisać istniejące pakiety FAR (N302, N328–N342, QW‑2191) do tego języka operatorów i które identyfikacje uznać za „kanoniczne”, a których nie.

