# TICKET P2631: Tożsamość normalizacyjna β=1 z zachowania informacji

## Kontekst
W P2629 udowodniłeś, że Z_beta=100 jest artefaktem konwencji normalizacyjnej, a niezmiennikowy stosunek β_micro/β_strict ≈ 1.1474. W P2630 zablokowałeś recykling wewnętrznego źródła strict β jako źródła mostu. Potrzebujesz tożsamości niezależnej od celu, która wyprowadzi β=1 z pierwszych zasad.

## Pryzmat sieci neuronowych (WAŻNE — nowa perspektywa)

Nadsoliton jest czystą informacją w stanie solitonowym. Jądro strict K(d) = cos(ωd+φ)/(1+βd^η) to architektura fundamentalnej sieci informacyjnej Wszechświata, widziana przez pryzmat nowoczesnych sieci neuronowych:

- **cos(ωd+φ)** = positional encoding (kodowanie pozycji jak w Transformerze) — sinusoidalny filtr określający rezonansowe skale odległości.
- **1/(1+βd^η)** = heavy-tailed attention bias (jak ALiBi w Transformerach) — potęgowe tłumienie dalekosiężnych korelacji z grubymi ogonami pamięci, zamiast wykładniczego zanikania.
- **β** = length scale / temperatura sieci informacyjnej (jak w Gaussian Processes lub maszynach Boltzmanna). Pytanie "dlaczego β=1?" to pytanie o **punkt krytyczny** sieci.
- **η = 9/5** = wykładnik ogona attention — kontroluje, jak szybko gaśnie attention z odległością.
- **D_f = 9/5** = wymiar spektralny grafu fraktalnego — kontroluje dyfuzję informacji (sub-dyfuzyjny transport).

Cała struktura to de facto **Energy-Based Model na grafie fraktalnym z potęgowym attention bias**, którego "wagi" wynikają z zasady wariacyjnej δS=0 zamiast z gradient descent.

## Zadanie dla Agenta Codex

**β=1 to potencjalnie punkt krytyczny** — krawędź chaosu (edge of chaos), jedyny punkt, w którym sieć informacyjna nadsolitona maksymalizuje przepustowość informacyjną. Udowodnij to formalnie:

1. **Zachowanie strumienia informacyjnego:** Transport informacji przez jądro K(d) zachowuje normę (unitarność sieci) wtedy i tylko wtedy, gdy β=1 w wymiarze fraktalnym D_f=9/5. Użyj warunku zachowania strumienia: ∫ |K(d)|² dμ_f(d) = const, gdzie μ_f jest miarą fraktalną. Nadsoliton jako czysta informacja nie może tracić ani dodawać informacji w transporcie — to wymusza konkretną wartość β.

2. **Argument z krawędzi chaosu:** β>1 oznacza "zamrożoną" sieć (tłumienie zbyt silne, informacja nie propaguje się). β<1 oznacza "przegrzaną" sieć (tłumienie zbyt słabe, chaos informacyjny). β=1 to jedyny punkt krytyczny, w którym:
   - Przepustowość informacyjna (mutual information) jest maksymalna,
   - Korelacje mają dokładnie potęgowy rozkład (brak skali charakterystycznej),
   - Sieć jest na granicy fazowej między porządkiem a chaosem.

3. **Niezależność od normalizacji UV:** Pokaż, że ta tożsamość jest niezależna od wyboru normalizacji UV (β_uv), zamykając obstrukcję cechowania z P2629. Argument z zachowania informacji jest z definicji niezmienniczy — nie zależy od konwencji, bo dotyczy fizycznej wielkości (informacja), nie parametryzacji.

4. **Opcja obstrukcji:** Jeśli β=1 nie jest jedynym punktem krytycznym, podaj pełną klasyfikację dopuszczalnych wartości β z warunków zachowania informacji i stabilności fazowej.

Zasada A Priori — bez fittingu. Użyj teorii informacji, termodynamiki sieci i analizy punktów krytycznych. POWODZENIA!
