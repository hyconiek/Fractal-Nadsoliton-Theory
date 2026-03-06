# A2 Supersoliton Matching

Status: `A2_EXECUTED_MINIMAL_MATCHING_BRANCH`
As of: `2026-03-06`

## Cel

Wykonac pierwszy rzeczywisty matching miedzy `A1` i tlem supersolitonowym, ale tylko na minimalnej galezi:
- `single-foundation`,
- `gauge-off`,
- `metric-spectator`.

To nie jest jeszcze pelny sektor fizyczny. To jest najciensza galaz, na ktorej da sie uczciwie zapisac rownania Eulera-Lagrange'a i rozdzielic:
- warunki wymuszone przez sam matching,
- warunki pozostawione swobodne,
- warunki zalezace od wyboru gauge / reprezentacji tla.

## Wykonana galaz A2

### Minimalny ansatz tla

Na wykonanej galezi przyjmujemy:

1. pole podstawowe:
   - `Psi^A(x) = rho(r) n^A(Omega)`,
   - gdzie `rho(r)` jest profilem radialnym,
   - `n^A(Omega)` niesie czesc orientacyjno-topologiczna i moze realizowac niezerowy ladunek topologiczny.
2. pole porzadku:
   - `Phi(x) = phi(r)`.
   - na poziomie ontologii programu `phi(r)` jest warstwa efektywnego porzadku wyekstrahowana z `Psi`, a nie drugim wspolfundamentalnym bytem.
3. sektor gauge:
   - `A_mu^I = 0` na galezi minimalnej,
   - sektor gauge pozostaje dopuszczony przez `A1`, ale nie jest jeszcze aktywowany w samym matchingu.
4. sektor metryczny:
   - `g_{mu nu} = eta_{mu nu}` na galezi minimalnej,
   - backreaction grawitacyjny zostaje odlozony poza `A2`.

## Zredukowana akcja radialna

Po podstawieniu minimalnego tla, matching redukuje sie do akcji:

```text
S_red = Omega_(d-1) int dt dr [
  r^(d-1) ( 1/2 K(rho,phi) (rho')^2 + 1/2 (phi')^2 - V_eff(rho,phi) )
  - r^(d-3) T_top(rho,phi; n, B)
]
```

Gdzie:
- `K(rho,phi)` jest efektywna metryka pola `G_AB` wzdluz tangensa tla,
- `V_eff(rho,phi) = V(rho,phi) + U_mix(rho,phi)`,
- `T_top` zbiera czesc katowo-topologiczna pochodzaca z `n^A(Omega)` i ewentualnego ladunku `B`.

To jest kluczowy punkt `A2`:
- matching nie wyznacza calego `G_AB`,
- matching nie wyznacza osobno `V` i `U_mix`,
- matching wyznacza tylko warunki na ich wartosci i pochodne wzdluz wykonanej galezi tla.

## Zredukowane rownania Eulera-Lagrange'a

Na galezi minimalnej dostajemy dwa jawne warunki matchingowe:

### Rownanie dla `rho`

```text
d/dr [ r^(d-1) K(rho,phi) rho' ]
- r^(d-1) [ 1/2 (partial_rho K) (rho')^2 - partial_rho V_eff ]
+ r^(d-3) partial_rho T_top
= 0
```

### Rownanie dla `phi`

```text
d/dr [ r^(d-1) phi' ]
- r^(d-1) [ 1/2 (partial_phi K) (rho')^2 - partial_phi V_eff ]
+ r^(d-3) partial_phi T_top
= 0
```

Interpretacja:
- `K` musi byc dodatnia i regularna wzdluz trajektorii tla,
- `V_eff` musi dopuszczac asymptotyczne prozniowe minimum,
- `T_top` nie moze rozwalac regularnosci rdzenia i musi gasnac dostatecznie szybko przy `r -> infinity`.

## Warunki brzegowe wymuszone przez matching

### Rdzen

Matching wymusza regularnosc w `r = 0`:
- `rho'(0) = 0`,
- `phi'(0) = 0`,
- `T_top` musi pozostawac skonczone,
- jesli wybrana galaz ma niezerowy ladunek topologiczny, to rdzen musi byc zgodny z degeneracja orientacyjna `n^A`.

### Nieskonczonosc

Finite-energy branch wymusza:
- `rho(r) -> rho_inf`,
- `phi(r) -> phi_inf`,
- `(rho_inf, phi_inf)` musi byc punktem krytycznym `V_eff`,
- `partial_rho V_eff(rho_inf, phi_inf) = 0`,
- `partial_phi V_eff(rho_inf, phi_inf) = 0`,
- `T_top / r^2 -> 0`.

## Tabela `forced / optional / gauge-choice-dependent`

| Obiekt | Klasyfikacja | Znaczenie |
|---|---|---|
| istnienie radialnych profili `rho(r), phi(r)` | `forced` | bez nich nie ma tla supersolitonowego |
| dodatni i regularny `K(rho,phi)` wzdluz tla | `forced` | inaczej matching traci sens jako baza do `A3` |
| istnienie asymptotycznego minimum `V_eff` | `forced` | potrzebne do finite-energy branch |
| zanik `T_top` na nieskonczonosci | `forced` | potrzebny do skonczonej energii |
| osobny rozklad `V` vs `U_mix` | `optional` | `A2` widzi tylko `V_eff = V + U_mix` wzdluz tla |
| pelna postac tensorowa `G_AB` poza trajektoria tla | `optional` | `A2` wymaga tylko projekcji na kierunek tla |
| aktywny profil `A_mu^I` | `optional` | minimalna galaz przyjmuje `A_mu^I = 0` |
| backreaction `g_{mu nu}` | `optional` | odlozony do pozniejszego bridge |
| orientacja `n^A` w przestrzeni wewnetrznej | `gauge-choice-dependent` | mozna ja ustalac konwencja lub gauge choice |
| `A_mu^I = 0` vs `A_mu^I` pure gauge | `gauge-choice-dependent` | obie realizacje sa zgodne z minimalna gala, dopoki nie niosa fizycznego fluxu |
| wybor wspolrzednych przy wlaczeniu metryki | `gauge-choice-dependent` | nie nalezy mylic z nowa trescia dynamiczna |

## Matching funkcjonalow z A1

### `G_AB`

- `forced`:
  - projekcja `G_AB` na tangent tla musi dawac dodatni `K(rho,phi)`,
  - brak osobliwosci na rdzeniu i na prozniowej granicy.
- `optional`:
  - cala poprzeczna geometria field-space pozostaje nieustalona.
- `gauge-choice-dependent`:
  - baza wewnetrzna, w ktorej zapisujemy `n^A`, nie jest fizycznie wyprowadzona przez samo `A2`.

### `V` i `U_mix`

- `forced`:
  - ich suma `V_eff` musi dopuszczac punkt krytyczny zgodny z warunkami brzegowymi.
- `optional`:
  - sam matching nie rozdziela jednoznacznie, co jest `V`, a co `U_mix`.

### `Z_IJ`

- `forced`:
  - brak wymuszenia na galezi `gauge-off`, poza tym ze ewentualna przyszla aktywacja nie moze psuc regularnosci.
- `optional`:
  - konkretna postac `Z_IJ` nie jest jeszcze wyznaczona.
- `gauge-choice-dependent`:
  - przejscie z `A_mu^I = 0` do `A_mu^I` pure gauge nie niesie jeszcze nowej tresci dynamicznej.

### `M_eff` i `Lambda_eff`

- `forced`:
  - na galezi minimalnej brak dodatkowego wymuszenia poza niesprzecznoscia z przyszlym wlaczeniem backreaction.
- `optional`:
  - konkretna postac pozostaje otwarta.

## Co A2 rzeczywiscie zamyka

`A2` zamyka tylko tyle:
- istnieje spójna minimalna galaz matchingu dla tla `single-foundation supersoliton`,
- wiadomo, jakie rownania i warunki brzegowe trzeba przeniesc do `A3`,
- wiadomo, ktore funkcjonaly sa rzeczywiscie widziane przez matching, a ktore pozostaja jeszcze swobodne.

## Co pozostaje otwarte po A2

- czy dla konkretnego wyboru `K`, `V_eff` i `T_top` rzeczywiscie istnieje rozwiazanie globalne,
- czy galaz `gauge-on` jest potrzebna do stabilizacji,
- czy topologiczna czesc `n^A(Omega)` daje jedyny dopuszczalny ladunek,
- jak wyglada pelny operator drugiej wariacji,
- czy po wlaczeniu fermionow i grawitacji matching pozostaje zgodny.

## Anti-overclaim

`A2` nie twierdzi, ze sam matching daje juz:
- spinory,
- konkretna grupe gauge,
- globalna stabilnosc,
- RG closure,
- GR-limit,
- finalny lagranzian SM+GR.

## Produkt etapu

- jawny minimalny ansatz tla,
- jawna zredukowana akcja radialna,
- podstawione rownania E-L dla `rho` i `phi`,
- tabela `forced / optional / gauge-choice-dependent`,
- lista wejsc do `A3`.

## Nastepny krok

Naturalnym nastepnym ruchem jest `A3`:
- policzyc druga wariacje wokol wykonanej galezi tla,
- rozdzielic `zero / gauge / physical modes`,
- zbudowac operator fluktuacji w dokladnie tej samej galezi, ktora zostala wykonana w `A2`.
