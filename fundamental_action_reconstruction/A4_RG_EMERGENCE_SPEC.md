# A4 RG Emergence Spec

Status: `A4_EXECUTED_ONE_STEP_MINIMAL_BRANCH_COARSE_GRAINING`
As of: `2026-03-07`

## Cel

Sprawdzic, czy z operatora fluktuacji z `A3` i jawnego rozdzialu skal wynika pierwszy uczciwy krok Wilsonowskiego coarse-graining, a nie tylko recznie dosztukowany proxy-RG.

Etap `A4` jest wykonany tylko na tej samej galezi minimalnej:
- `single-foundation`,
- `gauge-off`,
- `metric-spectator`.

## Ontologiczna wskazowka etapu

Ten etap jest prowadzony przy tej samej konstrukcyjnej ontologii:
- `Psi` jest jedynym polem traktowanym jako fundamentalne,
- `Phi` jest efektywna warstwa porzadku zwiazana z `Psi`,
- sektor gauge i metryczny pozostaja warstwami emergentnymi albo efektywnymi,
- to nadal nie jest theorem-level closure ontologii jednego nadsolitonu.

## Canonical fractal scaling supplement

Jesli `A4` ma byc RG-emergence layer dla tej samej ontologii FIN, to shell
integration nie moze juz byc opisywana tak, jakby podloze bylo zwyklym gladkim
substratem bez jawnego parametru fraktalnego.

Na obecnym etapie wolno wpisac tylko canonical-ontology-supported supplement:

```text
D_f ≡ alpha_geo ≡ 4 ln 2 ≈ 2.7726
beta_tors ≈ 0.01
```

Interpretacja:
- `D_f` nosi substratowy ciezar skalowania dla coarse-graining,
- `beta_tors` niesie warstwe tlumienia miedzy-layer / torsion-damping,
- nadal nie jest to globalny RG theorem ani strict derivation shell measure.

## Kernel source classification after `K1/K2/F2`

Po `K1`, `K2`, `P47/N50` i `F2` `A4` musi odroznic dwie role kernela:

1. ontologiczna warstwa action-first:
   `D_f / alpha_geo / beta_tors`,
2. pozniejszy strict working kernel:
   `K_strict_gate`.

Na current repo state wolno:
- sledzic shell scaling przez canonical-ontology-supported warstwe
  `D_f / beta_tors`,
- uzywac `K_strict_gate` tylko jako downstream operational control,
  benchmark albo consistency target.

Na current repo state nie wolno:
- traktowac `QW-2049` kernela jako dowodu, ze shell measure fraktalna jest juz
  wyprowadzona,
- milczaco podstawic `K_strict_gate` w miejsce ontologicznej warstwy
  `D_f / alpha_geo / beta_tors`.

## Wejscie z A3

`A3` dostarczylo:
- operator fizyczny `O_phys = - d/dr [ K_2(r) d/dr ] + M_2(r)`,
- split `zero / gauge / physical modes`,
- warunek, ze przyszle claimy o stabilnosci i RG wolno stawiac dopiero po projekcji na kanal fizyczny.

W `A4` pracujemy tylko na fizycznym podprzestrzennym pakiecie modow.

## Jednokrokowy coarse-graining

Rozdzielamy fizyczne fluktuacje:

```text
xi = xi_< + xi_>
```

gdzie:
- `xi_<` zbiera mody o skali `|p| < mu / b`,
- `xi_>` zbiera mody z shellu `mu / b <= |p| <= mu`,
- `b > 1` jest jednym krokiem coarse-graining.

Po scalkowaniu shellu `xi_>` zapisujemy:

```text
S_eff[xi_<] =
  S_phys[xi_<]
  + 1/2 Tr_shell^(D_f,beta_tors) log O_phys
  + Delta S_local
  + Delta S_EFT
```

Interpretacja:
- `S_phys` to wejscie z minimalnego kernela `A3`,
- `Tr_shell^(D_f,beta_tors) log O_phys` jest pierwszym uczciwym krokiem Wilsonowskiego przeplywu z jawnie zaznaczona warstwa fraktalnego skalowania i tlumienia,
- `Delta S_local` zbiera lokalne renormalizacje operatorow juz obecnych,
- `Delta S_EFT` zbiera wyzsze operatory generowane przez coarse-graining.

Na tym etapie wolno jeszcze tylko symbolicznie wskazac wage shellu:

\[
d\mu_{\mathrm{shell}}^{(f)} \sim p^{D_f-1} \, dp \, d\Omega_f
\cdot w_{\mathrm{tors}}(\beta_{\mathrm{tors}}),
\]

bez claimu, ze repo policzylo juz unikalna mikro-derewacje tej miary albo
globalny wieloskalowy przeplyw.

## Co jest emergentne, co wstawione recznie, a co nadal otwarte

| Obiekt | Status po A4 | Znaczenie |
|---|---|---|
| `K_tan(mu)` | `emergent` | efektywna tangensowa metryka kernela po shell integration |
| `H_V(mu)` | `emergent` | efektywny Hessian potencjalu na wykonanej galezi |
| `C_top(mu)` | `emergent` | efektywny topologiczny skladnik radialny po coarse-graining |
| `c_4(mu), c_6(mu), ...` | `emergent` | wyzsze lokalne operatory generowane jako ogon EFT |
| predeklarowane `Delta L_EFT` | `inserted_by_hand` | pozostaje jawnie wpisanym miejscem na wyzsze operatory |
| `Z_IJ(mu)` | `unresolved` | sektor gauge nie jest aktywny na tej galezi |
| `M_eff(mu), Lambda_eff(mu)` | `unresolved_but_fractal_scaling_relevant` | sektor grawitacyjny pozostaje spectator-only, ale nie wolno juz ignorowac `D_f/beta_tors` jako structural RG data |
| fermionowe sprzezenia biezace | `unresolved` | brak aktywnej galezi fermionowej w `A4` |

## Symboliczne beta-functions

Na poziomie wykonanej galezi wolno zapisac tylko symboliczne relacje:

```text
beta_K = Delta_K[ Tr_shell^(D_f,beta_tors) log O_phys ]
beta_H = Delta_H[ Tr_shell^(D_f,beta_tors) log O_phys ]
beta_top = Delta_top[ Tr_shell^(D_f,beta_tors) log O_phys ]
beta_c_n = canonical_part + shell_induced_part
```

To znaczy:
- `A4` pokazuje, skad ma sie brac running,
- running jest teraz jawnie obwarowany przez canonical-ontology-supported
  warstwe `D_f / beta_tors`,
- ale nie twierdzi jeszcze, ze policzono globalny RG package,
- nie twierdzi tez, ze running jest juz unikalny poza wykonana minimalna gala.

## Relevant / marginal / irrelevant

Na wykonanej galezi dostajemy minimalna klasyfikacje:

| Klasa | Kandydaci |
|---|---|
| `relevant` | masowe przesuniecia efektywnego Hessianu, skladniki prozniowe, niskowymiarowe deformacje radialne |
| `marginal` | dwupochodne skladniki bazowe w kanale fizycznym, wybrane sprzezenia bezwymiarowe jesli przetrwaja aktywacje dalszych sektorow |
| `irrelevant` | operatory wyzszego rzedu `c_4, c_6, ...` tlumione przez skale odciecia |

Ta tabela jest na razie lokalna dla wykonanej galezi i nie moze byc sprzedawana jako finalna klasyfikacja calej teorii.

## Co A4 rzeczywiscie zamyka

`A4` zamyka tylko tyle:
- istnieje jawny jednokrokowy schemat Wilsonowskiego coarse-graining dla minimalnego kernela z `A3`,
- wiadomo, ktore elementy runningu sa rzeczywiscie emergentne na tej galezi,
- wiadomo, ktore elementy pozostaja nadal wstawione recznie albo nierozstrzygniete,
- powstaje poprawny punkt wejscia do dalszego mostu `A5`.

## Co pozostaje otwarte po A4

- globalny, wieloskalowy RG closure,
- aktywacja sektora gauge w runningu,
- aktywacja sektora fermionowego,
- running sektora grawitacyjnego,
- powiazanie z historycznym `L12` bez dodatkowych theorem-level krokow,
- unikalnosc przeplywu poza wykonana minimalna gala.

## Anti-overclaim

`A4` nie twierdzi, ze:
- zamknieto globalny nonperturbative RG,
- z automatu zamknieto `L12`,
- z runningu juz wynika `SU(3)xSU(2)xU(1)`,
- zamknieto fermionowy albo grawitacyjny sektor przeplywu,
- uzyskano theorem-level/full-closure PASS.

Brak kontrprzykladu albo lokalna zgodnosc shell integration nie sa traktowane jako dowod globalnego closure.

## Produkt etapu

- jawny jednokrokowy schemat coarse-graining,
- efektywna akcja po jednym kroku RG,
- symboliczne beta-functions dla `K_tan`, `H_V`, `C_top`, `c_n`,
- tabela `emergent / inserted by hand / unresolved`,
- lista granic, ktorych nie wolno nadinterpretowac.

## Nastepny krok

Naturalnym kolejnym ruchem jest `A5`:
- rozdzielic droge spinor-emergent od minimal spin-bundle extension,
- ustalic, ktore dane z `A1..A4` rzeczywiscie wspieraja sektor fermionowy,
- i dopiero potem wracac do mostu SM+GR.
