# Fundamental Action Reconstruction

Status: `PROGRAM_OPENED_A1_COMPLETED_NO_FULL_CLOSURE_CLAIM`
As of: `2026-03-06`

Ten katalog otwiera rownolegly tor konstrukcyjny poza drabinka `L5/L12`.

Cel:
- zbudowac kandydacki lagranzian od strony dzialania,
- wymusic zgodnosc z tlem supersolitonowym,
- przejsc do analizy kernela fluktuacji,
- sprawdzic, czy RG wynika naturalnie z tej konstrukcji,
- dopiero potem wracac do spinorow, cechowania i mostu SM+GR.

Ten tor nie zastepuje drabinki. Drabinka zostaje jako:
- audit blockerow,
- anti-overclaim,
- walidacja lokalnych twierdzen provider-layer.

Ten tor jest glowna sciezka konstrukcyjna.

## Dlaczego ten tor jest lepszy od dalszego iterowania drabinki

Drabinka bardzo dobrze:
- izoluje aktualny blocker,
- pilnuje, gdzie proof chain sie urywa,
- zapobiega falszywym claimom theorem-level/full-closure.

Ale sama drabinka nie jest naturalnym mechanizmem konstrukcji pelnego lagranzianu. Obecnie glownie przesuwa nazwy nierozladowanych twierdzen. Dlatego najlepsza sciezka poza drabinka to:
- `action-first` jako sciezka glowna,
- `geometry/spectral` jako warstwa pomocnicza,
- drabinka jako audit i walidacja.

## Program A1..A10

1. `A1`: minimal action ansatz
   - zapisac najmniejszy lokalny ansatz dzialania zgodny z Lorentzem, topologia i celem supersolitonowym.
2. `A2`: supersoliton matching
   - wymusic, by rownania Eulera-Lagrange'a dopuszczaly zadane tlo supersolitonowe.
3. `A3`: kernel analysis
   - policzyc druga wariacje dzialania i rozdzielic zero modes, gauge modes i physical modes.
4. `A4`: RG emergence
   - zdefiniowac coarse-graining i sprawdzic, czy Wilsonowski flow wynika z integracji modow.
5. `A5`: spinor route split
   - rozdzielic dwie dopuszczalne drogi: emergent spinor sector vs minimal spin-bundle extension.
6. `A6`: gauge reconstruction
   - ustalic, kiedy i jak z konstrukcji wynika grupa cechowania oraz anomaly constraints.
7. `A7`: positivity/unitarity package
   - domknac poprawne kryteria dla sektora bozonowego, fermionowego i lorentzowskiego.
8. `A8`: gravity bridge
   - wyprowadzic most do sektora grawitacyjnego i warunkow GR-limit.
9. `A9`: SM+GR effective reduction
   - sprawdzic, jaki low-energy effective lagranzian rzeczywiscie redukuje sie z konstrukcji.
10. `A10`: calibration boundary and anti-overclaim audit
   - oddzielic realna derivation od calibration i wystawic jawny raport granic.

## Aktualny status

- `A1`: wykonane jako warstwa spec/ansatz + wygenerowany manifest.
- `A2`: gotowe do wykonania jako matching-spec.
- `A3`: gotowe do wykonania jako kernel-spec.
- `A4`: gotowe do wykonania jako RG-spec.
- `A5..A10`: backlog programu.

## Twarde ograniczenia rygoru

- brak claimu `theorem-level PASS`,
- brak claimu `full-closure PASS`,
- brak claimu, ze `SU(3)xSU(2)xU(1)` zostalo juz wyprowadzone,
- brak claimu, ze fermiony Diraca lub GR sa juz domkniete z tego toru,
- brak claimu, ze kernel jest po prostu "dodatnio okreslony" bez rozroznienia sektorow.

## Prawidlowe warunki fizyczne do dalszych etapow

- bosonic Euclidean sector: coercivity / positivity po gauge-fixingu i po odjeciu zerowych modow,
- fermionic sector: self-adjointness / ellipticity / positivity of `D^dagger D` / reflection positivity,
- Lorentzian sector: well-posedness, hyperbolicity, bounded-below Hamiltonian, brak ghostow.

## Pliki w tym katalogu

- `A1_MINIMAL_ACTION_ANSATZ.md`
- `A2_SUPERSOLITON_MATCHING_SPEC.md`
- `A3_KERNEL_ANALYSIS_SPEC.md`
- `A4_RG_EMERGENCE_SPEC.md`
- `A5_SM_GR_BRIDGE_SPEC.md`
- `a1_minimal_action_ansatz.py`
- `generated/a1_minimal_action_ansatz_summary.json`
- `manifest_action_reconstruction.json`
