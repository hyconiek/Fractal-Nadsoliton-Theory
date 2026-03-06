# Fundamental Action Reconstruction

Status: `PROGRAM_PHASE1_COMPLETE_C16_MINIMAL_COEFFICIENT_CLASS_TABLE_REDUCED`
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

## Ontologiczna wskazowka programu

Program jest prowadzony pod robocza ontologia:
- warstwa pierwotna ma charakter informacyjny,
- jednym fundamentalnym obiektem konstrukcyjnym jest nadsoliton,
- `Phi`, sektor gauge i sektor metryczny sa na tym etapie traktowane jako warstwy efektywne albo emergentne.

To nie jest theorem-level closure ontologii jednego nadsolitonu.
To jest konstrukcyjny kierunek pracy dla toru `A1..A10`.

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
   - zintegrowac to, co jest juz realnie zamkniete lokalnie, i jawnie wyodrebnic terminalne blockery globalne.
8. `A8`: gravity bridge
   - zintegrowac tylko effective/scope-closed gravity layers i jawnie oddzielic fundament od scope.
9. `A9`: SM+GR effective reduction
   - zlozyc wykonane warstwy material/gauge/gravity do jednej uczciwej warstwy effective.
10. `A10`: calibration boundary and anti-overclaim audit
   - finalnie oddzielic realna derivation od calibration i wystawic jawny raport granic.

## Drugi cykl

1. `B1`: mode uniqueness minimal extra structure audit
   - zawezic blocker `QW-2191` do jednego pytania fizycznego o pochodzenie selektora.
2. `B2`: internal orientation datum source audit
   - sprawdzic, czy strict core rzeczywiscie zawiera juz wewnetrzne zrodlo selectora, czy tylko control-route.
3. `B3`: topological selector bridge packet
   - zbudowac minimalny packet derivation `topological sign / FR branch -> selector`, bez udawania discharge.
4. `B4`: minimal sigma_int candidate
   - wskazac jeden kanoniczny kandydat `sigma_int`, bez udawania strict derivation.
5. `B5`: sigma_int local stability audit
   - sprawdzic, czy kandydat ma choc lokalna stabilnosc deformacyjna i jak daleko jest od full gauge-safety.
6. `B6`: sigma to selector factorized bridge
   - sprawdzic, czy `sigma_int_candidate` da sie uczciwie osadzic jako residualny `Z2` datum w selector route bez falszywego discharge.
7. `B7`: factorized selector mode scaffold compatibility audit
   - sprawdzic, czy factorized bridge jest zgodny z `QW-2190`, `QW-2191` i granicami `A6`.
8. `B8`: selector track anti-overclaim audit
   - zamknac mini-pakiet `B3_O1..B3_O5` bez falszywego PASS i z jawna lista residualnych blockerow.
9. `C1`: narrow foundational blocker selection
   - wskazac jeden dominujacy waski blocker po selector-track, zamiast trzymac zbyt szerokie `uniqueness open`.
10. `C2`: selector family internal origin packet
   - sprawdzic, czy pochodzenie rodziny `J_ab` daje sie przynajmniej zredukowac warunkowo do mniejszych obligacji.
11. `C3`: internal reference pair candidate from mode scaffold
   - sprawdzic, czy `QW-2190` zawiera juz techniczny kandydat pary referencyjnej, nawet jesli nie jest jeszcze fizycznym orientation datum.
12. `C4`: local quadratic mismatch kinematic reduction
   - sprawdzic, czy forma `J_ab` wynika juz kinematycznie z orbity rotacyjnej i lokalnego kosztu kwadratowego, nawet jesli metryka tego kosztu nie jest jeszcze wyprowadzona.
13. `C5`: projected Hessian selector-metric bridge
   - sprawdzic, czy selector family jest orbitalna forma projected second variation, nawet jesli ta projekcja nie jest jeszcze explicite wyeksportowana.
14. `C6`: projected second variation source audit
   - sprawdzic, czy strict core zawiera juz packet-ready komponenty projekcji drugiej wariacji na kandydacka plaszczyzne orientacji, nawet jesli nie zawiera jeszcze jawnych eksportow.
15. `C7`: mode pair to orientation slice schema packet
   - sprawdzic, czy istnieje juz class-level schema slownika `mode pair -> orientation slice`, nawet jesli brak jeszcze basis-level eksportu.
16. `C8`: projected block positivity descent audit
   - sprawdzic, czy dodatniosc projected block moze schodzic z juz certyfikowanego host-operatora, nawet jesli brak jeszcze relacji kompresji.
17. `C9`: action-origin host carrier audit
   - sprawdzic, czy host-operator z `QW-2186` i orientation slice z `C7` maja juz wspolny action-origin carrier, nawet jesli brak jeszcze jawnej identyfikacji host-to-Hessian i restrykcji do slice.
18. `C10`: Psi-sector host identification audit
   - sprawdzic, czy strict core ma juz packet-ready schema `QW-2186 host -> canonical Psi-sector quadratic carrier`, nawet jesli brak jeszcze konkretnego block-level matchingu.
19. `C11`: Psi-sector block extraction audit
   - sprawdzic, czy strict core ma juz packet-ready schema konkretnego `Psi-sector quadratic block`, nawet jesli brak jeszcze jawnego extraction/export package.
20. `C12`: minimal Psi-block extraction packet
   - sprawdzic, czy strict core ma juz reprezentatywny seed i jawny packet ekstrakcyjny dla `Psi-sector block`, nawet jesli brak jeszcze assembled submatrix.
21. `C13`: mode-basis control index-set audit
   - sprawdzic, czy strict core ma juz deterministyczny control index-set w bazie modowej, nawet jesli brak jeszcze transportu do canonical `Psi` basis.
22. `C14`: control mode-to-Psi transport schema
   - sprawdzic, czy strict core ma juz jawny control transport `mode basis -> Psi basis`, nawet jesli brak jeszcze fizycznej kanonizacji tego mostu.
23. `C15`: control-only pullback submatrix packet
   - sprawdzic, czy strict core ma juz formalny packet assembly `M_control = T_control^T H_PsiPsi T_control`, nawet jesli brak jeszcze coefficient-filled canonical `H_PsiPsi`.
24. `C16`: minimal Psi-Hessian coefficient-class table
   - sprawdzic, czy strict core ma juz reprezentatywne rows coefficient-class dla canonical `H_PsiPsi`, nawet jesli brak jeszcze exhaustive `12 x 12` coefficient table.

## Aktualny status

- `A1`: wykonane jako warstwa spec/ansatz + wygenerowany manifest.
- `A2`: wykonane na minimalnej galezi matchingu `single-foundation / gauge-off / metric-spectator` + wygenerowany summary JSON.
- `A3`: wykonane jako operator drugiej wariacji + split `zero / gauge / physical modes` na tej samej galezi minimalnej.
- `A4`: wykonane jako jednokrokowy coarse-graining / RG-emergence layer na tej samej galezi minimalnej.
- `A5`: wykonane jako spinor route split z lokalnym audytem metodologicznym poprzednich badan w repo.
- `A6`: wykonane jako strict-core gauge reconstruction layer z jawnym blockerem unikalnosci.
- `A7`: wykonane jako strict-scope positivity/unitarity package z jawnymi terminalnymi obligacjami `L5_O1a_O1` i `L5_O1b_O1`.
- `A8`: wykonane jako strict-scope partial gravity bridge z jawnymi foundational blockers dla `G`, EH i full SM+GR reduction.
- `A9`: wykonane jako strict-scope partial SM+GR effective reduction.
- `A10`: wykonane jako finalny audit pierwszego cyklu programu.
- `B1`: wykonane jako pierwszy krok drugiego cyklu; blocker gauge uniqueness zawężony do pytania o internal selector.
- `B2`: wykonane jako drugi krok drugiego cyklu; w obecnym strict core nie znaleziono wyprowadzonego `internal orientation datum`.
- `B3`: wykonane jako trzeci krok drugiego cyklu; istnieje juz packet `B3_O1..B3_O5`, ale derivation pozostaje pending.
- `B4`: wykonane jako czwarty krok drugiego cyklu; istnieje juz jeden kanoniczny kandydat `sigma_int`, ale strict bridge pozostaje open.
- `B5`: wykonane jako piaty krok drugiego cyklu; lokalna stabilnosc kandydata jest wsparta, ale pelny gauge quotient pozostaje open.
- `B6`: wykonane jako szosty krok drugiego cyklu; istnieje pierwszy jawny factorized bridge `(sigma_int_candidate, J_ab family) -> theta*=0`, ale `sigma_int` alone nadal nie wyprowadza selectora.
- `B7`: wykonane jako siodmy krok drugiego cyklu; factorized bridge jest zgodny ze scaffoldem `QW-2190` i z obstrukcja `QW-2191`, ale tylko jako control-route overlay.
- `B8`: wykonane jako osmy krok drugiego cyklu; selector-track ma juz audit `no false pass`, a residualne blockery pozostaja jawne.
- `C1`: wykonane jako pierwszy krok trzeciego mikrocyklu; dominujacy waski blocker to brak internal derivation rodziny selectorow `J_ab`.
- `C2`: wykonane jako drugi krok trzeciego mikrocyklu; pochodzenie `J_ab` zostalo zredukowane warunkowo do dwoch mniejszych sub-blockerow: brak wewnetrznej pary referencyjnej i brak dodatniego lokalnego kosztu kwadratowego.
- `C3`: wykonane jako trzeci krok trzeciego mikrocyklu; `QW-2190` dostarcza techniczny kandydat pary referencyjnej `(c1,s1)` / `(c2,s2)`, ale jego fizyczne podniesienie pozostaje open.
- `C4`: wykonane jako czwarty krok trzeciego mikrocyklu; forma `J_ab(theta)=2(a+b)(1-cos theta)` zostaje zredukowana kinematycznie na orbicie `O(2)`, ale fizyczna identyfikacja dodatniej lokalnej metryki mismatch pozostaje open.
- `C5`: wykonane jako piaty krok trzeciego mikrocyklu; selector family zostaje zwiazana warunkowo z projected Hessianem na kandydackiej plaszczyznie orientacji, ale brak jeszcze explicite wycietej projekcji i jej certyfikatu dodatniosci.
- `C6`: wykonane jako szosty krok trzeciego mikrocyklu; strict core zawiera juz packet-ready komponenty dla projected second variation, ale nie zawiera jeszcze ani jawnej mapy `mode plane -> fluctuation subspace`, ani plane-specific positivity certificate.
- `C7`: wykonane jako siodmy krok trzeciego mikrocyklu; class-level schema slownika `mode pair -> orientation-related slice` jest juz jawna, ale brak basis-level eksportu pozostaje.
- `C8`: wykonane jako osmy krok trzeciego mikrocyklu; dodatniosc projected block zostaje zawężona do problemu jawnej relacji kompresji do host-operatora z certyfikatem `QW-2186`.
- `C9`: wykonane jako dziewiaty krok trzeciego mikrocyklu; compression blocker zostaje zawężony dalej do dwoch brakujacych eksportow: host-operator -> Psi-sector quadratic carrier oraz carrier -> orientation slice.
- `C10`: wykonane jako dziesiaty krok trzeciego mikrocyklu; host-identification blocker zostaje zawężony dalej do braku coefficient-level lub block-level matchingu z konkretnym Psi-sector quadratic Hessian blockiem.
- `C11`: wykonane jako jedenasty krok trzeciego mikrocyklu; block-matching blocker zostaje zawężony dalej do braku jawnego extraction/export package dla konkretnego `Psi-sector quadratic block`.
- `C12`: wykonane jako dwunasty krok trzeciego mikrocyklu; extraction/export blocker zostaje zawężony dalej do braku assembled `Psi x Psi` submatrix i coefficient table dla jawnie wybranego index-set.
- `C13`: wykonane jako trzynasty krok trzeciego mikrocyklu; index-set blocker zostaje zawężony dalej do braku transportu z control mode basis do canonical `Psi` basis oraz braku assembled submatrix po tym transporcie.
- `C14`: wykonane jako czternasty krok trzeciego mikrocyklu; transport blocker zostaje zawężony dalej do braku fizycznej kanonizacji control transportu i braku assembled submatrix po jego przyjeciu.
- `C15`: wykonane jako pietnasty krok trzeciego mikrocyklu; assembled-submatrix blocker zostaje zawężony dalej do braku coefficient-filled canonical `H_PsiPsi` dla control-only pullback oraz braku restriction do orientation slice.
- `C16`: wykonane jako szesnasty krok trzeciego mikrocyklu; coefficient-filling blocker zostaje zawężony dalej do braku exhaustive `12 x 12` coefficient table oraz braku restriction do orientation slice.

## Twarde ograniczenia rygoru

- brak claimu `theorem-level PASS`,
- brak claimu `full-closure PASS`,
- brak claimu, ze `SU(3)xSU(2)xU(1)` zostalo juz wyprowadzone,
- brak claimu, ze fermiony Diraca lub GR sa juz domkniete z tego toru,
- brak claimu, ze kernel jest po prostu "dodatnio okreslony" bez rozroznienia sektorow,
- brak claimu, ze `A7` domknelo globalne `L5`,
- brak claimu, ze `A8` domknelo foundational GR bridge,
- brak claimu, ze `A9` domknelo unified SM+GR theorem package,
- brak claimu, ze wykonanie `A10` oznacza full ToE closure,
- brak claimu, ze `B1` domknelo axiom-free uniqueness,
- brak claimu, ze `B2` znalazlo internal selector,
- brak claimu, ze `B3` rozladowalo topological-selector bridge,
- brak claimu, ze `B4` wyprowadzilo `sigma_int` theorem-level,
- brak claimu, ze `B5` rozladowalo `B3_O2`.
- brak claimu, ze `B6` rozladowalo `B3_O3`.
- brak claimu, ze `B7` rozladowalo `B3_O4`.
- brak claimu, ze `B8` oznacza uniqueness closure.
- brak claimu, ze `C1` rozladowalo blocker foundational.
- brak claimu, ze `C2` dalo internal derivation `J_ab`.
- brak claimu, ze `C3` rozladowalo `C2_B1`.
- brak claimu, ze `C4` rozladowalo `C2_B2`.
- brak claimu, ze `C5` znalazlo projected Hessian albo rozladowalo `C2_B2`.
- brak claimu, ze `C6` znalazlo eksport projekcji albo rozladowalo `C5_B1`.
- brak claimu, ze `C7` daje basis-level dictionary albo rozladowuje `C6_B1`.
- brak claimu, ze `C8` daje plane-specific positivity certificate albo rozladowuje `C6_B2`.
- brak claimu, ze `C9` identyfikuje juz `QW-2186` host z Psi-sector Hessianem albo rozladowuje `C8_B1`.
- brak claimu, ze `C10` znajduje juz konkretny Psi-sector block albo rozladowuje `C9_B1`.
- brak claimu, ze `C11` wydobywa juz konkretny block w postaci gotowej do matchingu albo rozladowuje `C10_B1`.
- brak claimu, ze `C12` ma juz assembled submatrix albo rozladowuje `C11_B1`.
- brak claimu, ze `C13` daje juz canonical `Psi` index-set albo rozladowuje `C12_B1`.
- brak claimu, ze `C14` daje juz fizycznie kanoniczny transport albo rozladowuje `C13_B1`.
- brak claimu, ze `C15` daje juz coefficient-filled `Psi x Psi` submatrix albo rozladowuje `C14_B2`.
- brak claimu, ze `C16` daje juz exhaustive `12 x 12` coefficient table albo rozladowuje `C15_B1`.

## Zasada korzystania z poprzednich badan

- strict-admissible internal references wolno traktowac jako material konstrukcyjny,
- legacy / exploratory studies wolno traktowac tylko jako heurystyke albo negatywna kontrole,
- obecny tor nie traktuje legacy korpusu jako dowodu theorem-level.
- axiom-augmented closures nie sa wlaczane do rdzenia strict bez jawnego oznaczenia.

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
- `A6_GAUGE_RECONSTRUCTION_SPEC.md`
- `A7_POSITIVITY_UNITARITY_PACKAGE_SPEC.md`
- `A8_GRAVITY_BRIDGE_SPEC.md`
- `A9_SM_GR_EFFECTIVE_REDUCTION_SPEC.md`
- `A10_CALIBRATION_BOUNDARY_AND_ANTI_OVERCLAIM_AUDIT.md`
- `B1_MODE_UNIQUENESS_MINIMAL_EXTRA_STRUCTURE_AUDIT.md`
- `B2_INTERNAL_ORIENTATION_DATUM_SOURCE_AUDIT.md`
- `B3_TOPOLOGICAL_SELECTOR_BRIDGE_PACKET.md`
- `B4_MINIMAL_SIGMA_INT_CANDIDATE.md`
- `B5_SIGMA_INT_LOCAL_STABILITY_AUDIT.md`
- `B6_SIGMA_TO_SELECTOR_FACTORIZED_BRIDGE.md`
- `B7_FACTORIZED_SELECTOR_MODE_SCAFFOLD_COMPATIBILITY_AUDIT.md`
- `B8_SELECTOR_TRACK_ANTI_OVERCLAIM_AUDIT.md`
- `C1_NARROW_FOUNDATIONAL_BLOCKER_SELECTION.md`
- `C2_SELECTOR_FAMILY_INTERNAL_ORIGIN_PACKET.md`
- `C3_INTERNAL_REFERENCE_PAIR_CANDIDATE_FROM_MODE_SCAFFOLD.md`
- `C4_LOCAL_QUADRATIC_MISMATCH_KINEMATIC_REDUCTION.md`
- `C5_PROJECTED_HESSIAN_SELECTOR_METRIC_BRIDGE.md`
- `C6_PROJECTED_SECOND_VARIATION_SOURCE_AUDIT.md`
- `C7_MODE_PAIR_TO_ORIENTATION_SLICE_SCHEMA_PACKET.md`
- `C8_PROJECTED_BLOCK_POSITIVITY_DESCENT_AUDIT.md`
- `C9_ACTION_ORIGIN_HOST_CARRIER_AUDIT.md`
- `C10_PSI_SECTOR_HOST_IDENTIFICATION_AUDIT.md`
- `C11_PSI_SECTOR_BLOCK_EXTRACTION_AUDIT.md`
- `C12_MINIMAL_PSI_BLOCK_EXTRACTION_PACKET.md`
- `C13_MODE_BASIS_CONTROL_INDEX_SET_AUDIT.md`
- `C14_CONTROL_MODE_TO_PSI_TRANSPORT_SCHEMA.md`
- `C15_CONTROL_ONLY_PULLBACK_SUBMATRIX_PACKET.md`
- `C16_MINIMAL_PSI_HESSIAN_COEFFICIENT_CLASS_TABLE.md`
- `a1_minimal_action_ansatz.py`
- `a2_supersoliton_matching.py`
- `a3_kernel_analysis.py`
- `a4_rg_emergence.py`
- `a5_spinor_route_split.py`
- `a6_gauge_reconstruction.py`
- `a7_positivity_unitarity_package.py`
- `a8_gravity_bridge.py`
- `a9_sm_gr_effective_reduction.py`
- `a10_calibration_boundary_and_anti_overclaim_audit.py`
- `b1_mode_uniqueness_minimal_extra_structure_audit.py`
- `b2_internal_orientation_datum_source_audit.py`
- `b3_topological_selector_bridge_packet.py`
- `b4_minimal_sigma_int_candidate.py`
- `b5_sigma_int_local_stability_audit.py`
- `b6_sigma_to_selector_factorized_bridge.py`
- `b7_factorized_selector_mode_scaffold_compatibility_audit.py`
- `b8_selector_track_anti_overclaim_audit.py`
- `c1_narrow_foundational_blocker_selection.py`
- `c2_selector_family_internal_origin_packet.py`
- `c3_internal_reference_pair_candidate_from_mode_scaffold.py`
- `c4_local_quadratic_mismatch_kinematic_reduction.py`
- `c5_projected_hessian_selector_metric_bridge.py`
- `c6_projected_second_variation_source_audit.py`
- `c7_mode_pair_to_orientation_slice_schema_packet.py`
- `c8_projected_block_positivity_descent_audit.py`
- `c9_action_origin_host_carrier_audit.py`
- `c10_psi_sector_host_identification_audit.py`
- `c11_psi_sector_block_extraction_audit.py`
- `c12_minimal_psi_block_extraction_packet.py`
- `c13_mode_basis_control_index_set_audit.py`
- `c14_control_mode_to_psi_transport_schema.py`
- `c15_control_only_pullback_submatrix_packet.py`
- `c16_minimal_psi_hessian_coefficient_class_table.py`
- `generated/a1_minimal_action_ansatz_summary.json`
- `generated/a2_supersoliton_matching_summary.json`
- `generated/a3_kernel_analysis_summary.json`
- `generated/a4_rg_emergence_summary.json`
- `generated/a5_spinor_route_split_summary.json`
- `generated/a6_gauge_reconstruction_summary.json`
- `generated/a7_positivity_unitarity_package_summary.json`
- `generated/a8_gravity_bridge_summary.json`
- `generated/a9_sm_gr_effective_reduction_summary.json`
- `generated/a10_calibration_boundary_and_anti_overclaim_audit_summary.json`
- `generated/b1_mode_uniqueness_minimal_extra_structure_audit_summary.json`
- `generated/b2_internal_orientation_datum_source_audit_summary.json`
- `generated/b3_topological_selector_bridge_packet_summary.json`
- `generated/b4_minimal_sigma_int_candidate_summary.json`
- `generated/b5_sigma_int_local_stability_audit_summary.json`
- `generated/b6_sigma_to_selector_factorized_bridge_summary.json`
- `generated/b7_factorized_selector_mode_scaffold_compatibility_audit_summary.json`
- `generated/b8_selector_track_anti_overclaim_audit_summary.json`
- `generated/c1_narrow_foundational_blocker_selection_summary.json`
- `generated/c2_selector_family_internal_origin_packet_summary.json`
- `generated/c3_internal_reference_pair_candidate_from_mode_scaffold_summary.json`
- `generated/c4_local_quadratic_mismatch_kinematic_reduction_summary.json`
- `generated/c5_projected_hessian_selector_metric_bridge_summary.json`
- `generated/c6_projected_second_variation_source_audit_summary.json`
- `generated/c7_mode_pair_to_orientation_slice_schema_packet_summary.json`
- `generated/c8_projected_block_positivity_descent_audit_summary.json`
- `generated/c9_action_origin_host_carrier_audit_summary.json`
- `generated/c10_psi_sector_host_identification_audit_summary.json`
- `generated/c11_psi_sector_block_extraction_audit_summary.json`
- `generated/c12_minimal_psi_block_extraction_packet_summary.json`
- `generated/c13_mode_basis_control_index_set_audit_summary.json`
- `generated/c14_control_mode_to_psi_transport_schema_summary.json`
- `manifest_action_reconstruction.json`
