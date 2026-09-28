from pathlib import Path
import json,csv,hashlib,math
ROOT=Path('/mnt/data/fin_physical_bridge_campaign_01')
def load(p): return json.load(open(ROOT/p,encoding='utf-8'))
def dump(p,o): (ROOT/p).write_text(json.dumps(o,indent=2,ensure_ascii=False,sort_keys=True),encoding='utf-8')
def write(p,s): (ROOT/p).write_text(s.strip()+"\n",encoding='utf-8')
rep=load('PHYS-002/REPRODUCTION.json'); fp=load('PHYS-003/FINGERPRINT.json'); mods=load('PHYS-004/COUNTERMODELS.json'); pr=load('PHYS-007/PREREGISTRATION.json'); pred=load('PHYS-007/PREDICTIONS.json')
A=rep['A7']; rho=rep['rho3']
# Common epistemic footer
nonc='''# NONCONCLUSIONS\n\n- Wynik nie wyprowadza fizycznej konieczności A7 ani parametrów jej kernela.\n- Nie wyprowadza SI czasu, przestrzeni, QM, GR, Modelu Standardowego ani ToE.\n- Programowana realizacja A7 testuje realizację zadanego modelu, nie pochodzenie A7 w przyrodzie.\n- Protokół 333 (`theta=2`, `Tprep=4`) pozostaje niezmieniony i nie jest testowany ani retunowany w tej kampanii.\n- Task 335 nie został otwarty. Task 337 pozostaje WIP.\n- `Gamma=B4` nie jest promowane ponad aktualny scoped intake: pozostają braki enclosure/LP/scalar-minimum i pełnego replay.\n'''
# PHYS001
r1={'execution_status':'COMPLETED','scientific_verdict':'EXACT','claim_class':'EXACT','accepted':True,'enumeration_tests':load('PHYS-001/enumeration_tests.json'),'key_result':'Finite-N FIN Gibbs core is exactly a generalized 12-state matrix Curie-Weiss model; equivalently a discrete vector-spin mean-field model with 12 allowed vectors in R^7. This does not mean seven physical spatial dimensions.'}
dump('PHYS-001/RESULTS.json',r1)
write('PHYS-001/DICTIONARY.md',r'''# PHYS-001 — exact Curie–Weiss dictionary

## Dokładna równoważność
Dla etykiet `sigma_a in Z_12` definiujemy wyłącznie jako słownik modelowy

`H(sigma)=-(J_c/(2N)) sum_{a,b} A7[sigma_a,sigma_b] - h sum_a 1[sigma_a=0]`.

Przy `beta_th=1/(k_B T)`, `g=beta_th J_c` i `vartheta=beta_th h=kappa/N` miara kanoniczna daje dla liczebności `n` dokładnie

`pi(n) ∝ N!/prod_j n_j! * exp[(g/(2N)) n^T A7 n + vartheta n_0]`.

Współczynnik `1/2` usuwa podwójne liczenie par w sumie po `a,b`; `1/N` jest skalowaniem mean-field zapewniającym energię ekstensywną. Czynnik multinomialny jest dokładną liczbą mikrostanów o danych liczebnościach.

Po usunięciu aktualizowanej kopii `i`, `m=n-e_i`. Ponieważ A7 ma stałą przekątną, składnik samodiagonalny nie zależy od proponowanej etykiety i znika w normalizacji:

`q_j(m)=softmax_j[(g/N)(A7 m)_j + vartheta 1[j=0]]`.

Stąd count-rate dla heat-bath wynosi `nu n_i q_j(m)`; `nu=1` w repo jest wyborem jednostki generator-time, a nie sekundą.

## Generalized Curie–Weiss, Potts i vector-spin
To jest dokładny macierzowy/generalized q-state Curie–Weiss. Nie jest automatycznie najprostszym Potts Hamiltonianem `-J delta_{sigma_a,sigma_b}`, bo A7 nie ma postaci jednego centered Potts projector. Ponieważ `A7=X7 X7^T`, można pisać `xi_j=X7[j,:]` i

`H=-(J_c/(2N)) |sum_a xi_{sigma_a}|^2 - h n_0`.

Jest to więc także dyskretny wektorowy model mean-field z 12 dozwolonymi wektorami wewnętrznymi w `R^7`. `R^7` jest przestrzenią cech/spinów tego modelu, nie wyprowadzoną 7-wymiarową przestrzenią fizyczną.

## Gauge skali
Transformacja `A -> c A`, `J_c -> J_c/c` pozostawia Hamiltonian parowy bez zmian. Zatem sama skala A nie jest identyfikowalna bez konwencji dla `J_c`/`g`; kształt widma po normalizacji jest osobną informacją.

## Granica roszczenia
Równoważność miary Gibbsa jest matematyczna (M). Wybór heat-bath jest dodatkową regułą kinetyczną. Utożsamienie `J_c,h,T,nu` z wielkościami konkretnego urządzenia wymaga osobnego mostu fizycznego/kalibracji; target Hamiltonian nie jest przez to energią elektryczną sprzętu.
''')
dump('PHYS-001/MAPPING.json',{'sigma':'12-state internal label','n_j':'occupation count','p_j':'n_j/N','A7':'dimensionless matrix interaction structure','g':'beta_th*J_c','vartheta':'beta_th*h=kappa/N','N':'number of copies/elements','nu':'attempt frequency; repository nu=1 defines time unit','claim_levels':{'Hamiltonian/Gibbs/count degeneracy':'M exact','hardware energy':'not identified','physical temperature':'not identified unless calibrated','R7 feature space':'internal feature space, not physical space'}})
write('PHYS-001/REPORT.md',f'''# PHYS-001 REPORT\n\n**execution_status:** COMPLETED  \n**scientific_verdict:** EXACT\n\nRównoważność z generalized 12-state matrix Curie–Weiss przechodzi dokładnie. Test agregacji mikrostanów do count states dał maksymalny błąd {max(r1['enumeration_tests']['1']['max_aggregated_probability_error'],r1['enumeration_tests']['2']['max_aggregated_probability_error']):.3e} dla N=1,2. Kluczowa granica: Gibbs nie wybiera kinetyki, a siedem cech A7 nie oznacza siedmiu wymiarów przestrzeni fizycznej.\n''')
write('PHYS-001/NONCONCLUSIONS.md',nonc)
write('PHYS-001/REPLAY.md','''# REPLAY\n`python ../run_smallN_campaign.py`\n\nSprawdź `enumeration_tests.json`; test nie używa N>2 dla tego taska.''')
write('PHYS-001/proof_enumerate.py','''from pathlib import Path\nimport runpy\nrunpy.run_path(str(Path(__file__).parents[1]/"run_smallN_campaign.py"),run_name="__main__")\n''')

# PHYS002
r2={'execution_status':'COMPLETED_WITH_SOURCE_GAP','scientific_verdict':'PASS_SCOPED_RECONSTRUCTED_CORE__ORIGINAL_FIN_SON_UNAVAILABLE','claim_class':'NUMERICAL_EVIDENCE+EXACT_STRUCTURE','accepted_for_downstream_smallN':True,'rho3':rho,'A7':A,'source_gap':'Historical FIN son/fin_core.py and seed_checks.py referenced by roadmap were not tracked at baseline commit and could not be byte-for-byte audited. A fail-closed small-N module was reconstructed from versioned formulas/sources instead.','barrier_status':'No promotion: V_d4 treated local; Gamma=B4 remains unpromoted under scoped intake.'}
dump('PHYS-002/RESULTS.json',r2)
write('PHYS-002/CODE_REVIEW.md',f'''# PHYS-002 — reproduction and seed-code audit\n\n## Verdict\n`PASS_SCOPED_RECONSTRUCTED_CORE__ORIGINAL_FIN_SON_UNAVAILABLE`. Repozytoryjne, śledzone źródła wystarczają do niezależnego odtworzenia małego rdzenia, ale dwa historyczne pliki `FIN son` nie są dostępne w bazowym commicie, więc nie twierdzę, że zostały naprawione bajt-w-bajt.\n\n## Fail-closed checks\nA7: rank={A['rank_tol_1e-10']}, trace={A['trace']:.12g}, centering residual={A['centering_max_abs']:.3e}, diagonal spread={A['diag_spread']:.3e}, min eigenvalue={A['min_eigenvalue']:.3e}. Wszystkie trzy N=2 generatory przechodzą row-sum/stationarity/detailed-balance bez symetryzowania błędnego wejścia.\n\nKlasyfikacja modów używa projektorów/orbit C12, a nie etykiety pierwszego wektora z degenerowanego multipletu. Przy N=2,g=0 wszystkie 12 sektorów są obecne; maksymalna eigenpair residual jest {rep['sector_g0_N2']['max_eigenpair_residual']:.3e}.\n\nDla N=3, heat-bath, `g=G_FROZEN` sektor k=4 daje `rho={rho['computed']:.17g}`, różnica od fixture {rho['abs_error']:.3e}.\n\n`G_FROZEN` zgadza się z R09 `g_bal` do ~2e-15 i ma pochodzenie operacyjne (balans kanałów bariery), ale nie jest stałą fizyczną. `V_d4` nie jest używane jako dowód globalny. Aktualny scoped intake nadal wymaga dokładnych enclosure wejść, LP i scalar minimum przed promocją `Gamma=B4`.\n''')
write('PHYS-002/NONCONCLUSIONS.md',nonc+'''\n- Brak byte-identical `FIN son` oznacza, że wynik dotyczy odtworzonego, wersjonowanego kontraktu, nie certyfikacji oryginalnych dwóch plików.\n''')
write('PHYS-002/REPLAY.md','''# REPLAY\n`python -m pytest -q test_safe_seed.py`\n`python safe_seed.py --N 2 --g 5.145228719489142`\n`python ../run_smallN_campaign.py`\n\nModuł nie uruchamia kosztownych obliczeń przy imporcie; CLI ogranicza N do 1..4.''')

# PHYS003
maxder=max(v['abs_error'] for N in fp['derivative_checks'].values() for v in N.values())
r3={'execution_status':'COMPLETED','scientific_verdict':'EXACT_SMALL_N_FINGERPRINT_PASS','claim_class':'EXACT+NUMERICAL_EVIDENCE','accepted':True,'active_ratios':A['active_ratios_to_k3'],'max_derivative_replay_abs_error':maxder,'finite_g_max_remainder_N3_g0p2':fp['finite_g_remainder']['3']['0.2']['max_abs']}
dump('PHYS-003/RESULTS.json',r3)
write('PHYS-003/PROOF.md',r'''# PHYS-003 — static spectral fingerprint

Definiuj `z_k(sigma)=exp(2 pi i k sigma/12)` oraz
`S_k=(1/N) E |sum_a z_k(sigma_a)|^2`.

Przy `g=0` etykiety są niezależne i jednostajnie rozłożone. Diagonalne składniki `a=b` dają N, a dla `a!=b` średnia charakteru niezerowego k znika. Zatem `S_k(N,0)=1` dla k=1..11.

Dla miary `pi_g ∝ exp[(g/(2N)) n^T A n]`,
`d<E F>/dg = Cov_g(F,n^T A n/(2N))`.
W `g=0`, po rozwinięciu sum parowych i użyciu ortogonalności znaków Z12, przeżywają tylko pary odpowiadające temu samemu sektorowi Fouriera. Daje to dokładnie

`dS_k/dg|_0 = ((N-1)/(12N)) Lambda_k`.

Dla FIN `Lambda_1=Lambda_2=0`, zaś aktywne są k=3,4,5,6. Po usunięciu skali trzy niezależne liczby kształtu wynoszą:
`Lambda4/Lambda3=1.12142406146`, `Lambda5/Lambda3=1.17191711554`, `Lambda6/Lambda3=1.19413370400`.

Dla N=2 i stałej przekątnej dokładny histogram różnicy etykiet ma postać
`P(d|i=0)=softmax_d[(g/2) A[0,d]]`, więc log-odds mierzą `(g/2)(A[0,d]-A[0,d'])` bez przybliżenia małego g.

Ważne: zerowy *slope* k=1,2 nie oznacza `S_1=S_2=1` dla skończonego g; nieliniowe poprawki mogą być niezerowe. N=1 jest kontrolą zerową: wszystkie S_k=1 dla dowolnego A,g.
''')
write('PHYS-003/REPORT.md',f'''# PHYS-003 REPORT\n\n**execution_status:** COMPLETED  \n**scientific_verdict:** EXACT_SMALL_N_FINGERPRINT_PASS\n\nCentralne różniczkowanie dokładnych sum N=2,3 odtwarza analityczne nachylenia z maksymalnym błędem {maxder:.3e}. Dla N=3 przy g=0.2 maksymalna reszta względem liniowej odpowiedzi wynosi {fp['finite_g_remainder']['3']['0.2']['max_abs']:.3e}; dla testu N=2 można całkowicie ominąć rozwinięcie małego g, używając dokładnego 12-bin pair histogram.\n\nFingerprint statyczny zależy od A7, a nie od wyboru heat-bath/Metropolis/Barker, o ile wszystkie mają tę samą miarę stacjonarną.\n''')
write('PHYS-003/FINITE_G_ERROR_BUDGET.md','''# Finite-g error budget\nDla pochodnych użyto dokładnej enumeracji count states N=2,3. Przy g=0.2 maksymalna odchyłka od liniowej odpowiedzi jest ~1.04e-3 (N=2) i ~9.47e-4 (N=3). Pierwszy proponowany test walidacyjny używa jednak dokładnego N=2 pair histogramu, więc jego błąd „finite-g approximation” wynosi 0 w obrębie modelu; pozostają kalibracja i statystyka.''')
write('PHYS-003/enumerator.py','''from pathlib import Path\nimport sys,json\nsys.path.insert(0,str(Path(__file__).parents[1]/"PHYS-002"))\nimport safe_seed as s\n_,_,_,_,A=s.build_rank7()\nfor N in (1,2,3):\n st=s.compositions(N); pi=s.count_pi(st,0.2,A)\n print(N,{k:s.static_S(st,pi,N,k) for k in range(1,7)})\n''')
write('PHYS-003/NONCONCLUSIONS.md',nonc)
write('PHYS-003/REPLAY.md','''# REPLAY\n`python enumerator.py`\n`python ../run_smallN_campaign.py`\nSprawdź `FINGERPRINT.json` / `derivative_checks`.''')

# PHYS004
# Get dynamic k4 values
rows=list(csv.DictReader(open(ROOT/'PHYS-004/CONTRAST_MATRIX.csv')))
dyn={}
for x in rows:
 if x['kind']=='dynamic_sector' and x['sector_k']=='4': dyn[x['model']]=float(x['normalized_rate'])
r4={'execution_status':'COMPLETED','scientific_verdict':'PASS_STATIC_SPECIFICITY__DYNAMIC_NOT_A7_ONLY','claim_class':'NUMERICAL_EVIDENCE+EXACT_EQUILIBRIUM_INVARIANCE','accepted':True,'time_normalization':'each kinetic scaled so its N=2,g=0,k=4 rate equals 1','N2_G_FROZEN_k4_rates':dyn,'diagnosis':'same A,g,pi but materially different relaxation rates => dynamic rates are kinetic-rule dependent, not A7-only'}
dump('PHYS-004/RESULTS.json',r4)
dump('PHYS-004/CALIBRATION_CONTRACT.json',{'matrix_scale':'all countermodels trace-matched to FIN trace before comparison','g_comparison':'same dimensionless g after trace match','no_g_over_geq_matching':True,'static_g_points':[0.2,1.0,2.0,3.0],'dynamic_point':s.G_FROZEN if False else 5.145228719489142,'dynamic_clock_normalization':'normalize each kinetic by its N=2 g=0 k4 relaxation rate; no post-hoc fit to FIN rates'})
write('PHYS-004/REPORT.md',f'''# PHYS-004 — countermodels and kinetic robustness\n\n**execution_status:** COMPLETED  \n**scientific_verdict:** PASS_STATIC_SPECIFICITY__DYNAMIC_NOT_A7_ONLY\n\nPrzed porównaniem zamrożono trace A i użyto tego samego g. Zbudowano full centered Potts, flat rank-7 P7, trzy perturbacje trace-preserving oraz osobny 2% leakage do k=1,2.\n\nStatyczny N=2 pair histogram przy g=3 odróżnia FIN m.in. od full Potts (TV≈0.096999), flat P7 (TV≈0.011787) i dwóch 10% perturbacji aktywnych wag (TV≈0.0128–0.0148). Najbliższe zaplanowane sondy 5%/2% mają TV≈0.00442/0.00563 i są znacznie trudniejsze.\n\nDynamika nie jest A7-only. Po jawnej normalizacji zegara `rate_k4(g=0)=1`, przy tym samym A7 i g=G_FROZEN sektor k=4 ma: heat-bath={dyn['heat_bath']:.6f}, Metropolis={dyn['metropolis']:.6f}, Barker={dyn['barker']:.6f}. Każda kinetyka zachowuje tę samą pi, ale relaksuje inaczej.\n\nWniosek: statyka może testować kształt A7 niezależnie od kinetyki; R_k wymaga osobnego, jawnego kontraktu generatora.\n''')
write('PHYS-004/NONCONCLUSIONS.md',nonc+'''\n- Nie istnieje na tej podstawie „uniwersalna dynamika FIN”; trzy dozwolone kinetyki są kontrprzykładem.\n''')
write('PHYS-004/REPLAY.md','''# REPLAY\n`python ../run_smallN_campaign.py`\nOdczytaj `CONTRAST_MATRIX.csv`; nie dopasowuj g po zobaczeniu wyników.''')

# PHYS005
r5={'execution_status':'COMPLETED','scientific_verdict':'NO_NEW_SOURCE','claim_class':'NO_GO/OPEN','accepted':True,'remaining_shape_parameters_after_D12_PSD_support_trace':3,'source_chain':'supplied strict kernel parameters -> circulant Laplacian spectrum -> manually selected support k=3,4,5,6 -> X7 -> A7','gaussian_mediator':'valid engineered parent representation, not independent source unless K and c_j are independently fixed/measured before A7 fingerprint','max_new_source_atoms':2}
dump('PHYS-005/RESULTS.json',r5)
dump('PHYS-005/SOURCE_PREMISES.json',{'nodes':[{'step':'strict kernel W(d)','status':'SUPPLIED','open':'physical origin of omega, phi, exponent/scale'},{'step':'Fourier/Laplacian spectrum','status':'DERIVED_FROM_SUPPLIED_KERNEL'},{'step':'support k=3,4,5,6','status':'SUPPLIED_SELECTION','open':'selection rule/necessity'},{'step':'four positive active eigenweights','status':'DERIVED_NUMERICALLY'},{'step':'trace normalization','status':'CONVENTION'},{'step':'three active spectral ratios','status':'FIN_SPECIFIC_SHAPE'},{'step':'A7=X7X7^T','status':'ALGEBRAIC_CONSTRUCTION'}],'candidate_source_atoms':[{'id':'MEDIATOR_MEASUREMENT','test':'independently measure/predict K and state couplings c_j, then predict A_ij=c_i^T K^-1 c_j without fitting fingerprint'},{'id':'SELECTION_RULE','test':'independent symmetry/physical rule predicts both k1,k2 cutoff and active eigenvalue ratios before FIN data'}]})
write('PHYS-005/SOURCE_GRAPH.md',r'''# PHYS-005 — A7 source and necessity ledger

`supplied strict kernel W(d)` → `circulant Laplacian Fourier spectrum` → **selected** support `{3,4,5,6}` → feature matrix `X7` → `A7=X7 X7^T`.

Pierwsza i pogrubiona strzałka nie są obecnie niezależnie fizycznie wyprowadzone. D12, PSD i rank 7 ograniczają klasę, ale po ustaleniu aktywnego supportu pozostają cztery dodatnie wagi, a po trace-normalization trzy niezależne proporcje. Zatem symmetry/rank/trace nie wybiera konkretnego tuple FIN.

Gaussian-parent jest poprawnym mostem konstrukcyjnym:
`H_parent=1/2 x^T K x - N^-1/2 x^T sum_a c_{sigma_a} - h n0`.
Eliminacja Gaussowskiego x przez completion of square daje pair matrix proporcjonalną do `c_i^T K^-1 c_j`. Jednak dowolna PSD A ma taką faktoryzację; wybór K,c po to, by odtworzyć A7, jest engineered realization. Staje się źródłem predykcyjnym dopiero wtedy, gdy K i c są ustalone z niezależnej fizyki/pomiarów przed oglądaniem A7 fingerprint.

**Verdict: NO_NEW_SOURCE.** Zachowujemy najwyżej dwa falsyfikowalne atomy: niezależny mediator measurement i niezależną selection rule. Bez nich nie dodajemy nowej „zasady fundamentalnej”.
''')
with open(ROOT/'PHYS-005/necessity_universality.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['feature','class','necessity_status']);w.writerows([['Gibbs/mean-field','standard statistical mechanics','not FIN-specific'],['q=12','FIN model choice','not independently necessary'],['support k=3,4,5,6','FIN-specific','source open'],['active eigenvalue ratios','FIN-specific','source open'],['Gaussian PSD mediator representation','generic for PSD matrices','engineered unless K,c independently sourced'],['metastability/FDT','standard mechanisms','not source evidence']])
write('PHYS-005/NONCONCLUSIONS.md',nonc+'''\n- Nie znaleziono niezależnego źródła A7; to terminalny `NO_NEW_SOURCE`, nie powód do wymyślenia nowego postulatu.\n''')
write('PHYS-005/REPLAY.md','''# REPLAY\n`python ../run_smallN_campaign.py` odtwarza widmo/trace/kontrmodele. Source verdict jest ledgerem zależności, nie wynikiem dopasowania.''')

# PHYS006
r6={'execution_status':'COMPLETED_DESIGN_ONLY','scientific_verdict':'PASS_SCOPED_REALIZATION_CONTRACT__NO_EXISTING_EXACT_Q12_DEVICE_CLAIM','claim_class':'DESIGN_ONLY','selected_platform':'programmable mixed-signal/electronic categorical sampler with explicit controller and physical stochastic source','hardware_execution_authorized':False,'laboratory_contact_authorized':False}
dump('PHYS-006/RESULTS.json',r6)
dump('PHYS-006/DEVICE_ASSUMPTIONS.json',{'required_operations':['12 valid categorical states per copy','maintain seven collective feature sums M','self-removal M-xi_old before logits','12 programmable logits ell_j','categorical draw calibrated to softmax','explicit asynchronous/serial update clock','microscopic label readout'],'not_demonstrated_by_current_literature':['one exact q=12 FIN node with arbitrary 12 logits','full arbitrary A7 leave-one-out generator','calibrated generator-time equivalence'],'Hamiltonian_status':'sampler target/objective unless a separate thermal-energy calibration is demonstrated','rng':'physical stochastic source desirable for physical realization; digital PRNG does not establish fundamental randomness'})
write('PHYS-006/DEVICE_MAPPING.md',r'''# PHYS-006 — physical realization contract

## Wybrana platforma
Najmniej dodatkowych założeń wymaga programowalny elektroniczny/mixed-signal sampler kategoryczny: kontroler przechowuje 12-state labels, liczy siedem pól zbiorowych `M=sum xi_sigma`, odejmuje własny `xi_old`, tworzy 12 logitów
`ell_j=(g/N) xi_j·(M-xi_old)+vartheta 1[j=0]`, a skalibrowany fizyczny element losowy realizuje kategorię softmax. Zegar aktualizacji jest jawny.

## Minimalny prototyp
N=1,g=0: kontrola jednostajnego 12-state losowania i bias/readout.  
N=2: jedna etykieta może być kotwicą; druga ma realizować dokładnie `P(j|i)=softmax[(g/2)A7[i,j]]`. To wystarcza do pierwszego statycznego pair testu bez bariery/metastability.

## Literatura a projekcja
FB-MOSFET Potts p-bits (Advanced Materials 2026, PMCID PMC12994325) demonstrują fizyczną stochastyczność, one-hot multi-state sampling i Boltzmann-like systemy, ale publikacja nie jest demonstracją dokładnego q=12 FIN/A7 leave-one-out sampler. Autorzy pokazują głównie q=4 hardware i kontrolowane wielostanowe probabilistyczne jednostki.  
https://pmc.ncbi.nlm.nih.gov/articles/PMC12994325/

Coupled CMOS ring-oscillator Potts machine (arXiv:2504.11376) jest przede wszystkim architekturą optimization/graph-coloring. Samo znajdowanie niskiej energii nie certyfikuje Gibbs distribution ani wymaganych conditional rates.  
https://arxiv.org/abs/2504.11376

## Kalibracja konieczna przed jakimkolwiek EV
Trzeba zmierzyć `P(new label | old environment)` dla reprezentatywnych 12-logitowych wektorów, całkowity attempt rate, korelacje RNG, latency/jitter i rzeczywisty schedule. Target H pozostaje funkcjonałem samplera, dopóki osobno nie wykazano fizycznej relacji `beta_th H` do energii/temperatury urządzenia.
''')
write('PHYS-006/CALIBRATION_PLAN.md','''# Calibration plan\n1. N=1,g=0: uniformity, readout-confusion matrix, RNG autocorrelation.\n2. Clamp one N=2 label; sweep pre-registered logit patterns and compare all 12 conditional probabilities to softmax.\n3. Verify self-removal and label permutation covariance.\n4. Measure attempt-rate distribution, delays and serial/parallel scheduling.\n5. Repeat on independent calibration records; validation records remain unopened.\n6. Go/no-go for static test: calibrated conditional prediction must fit within pre-registered TV envelope. Dynamic FIN test requires a stronger rate/schedule contract and is not authorized here.''')
with open(ROOT/'PHYS-006/FEASIBILITY.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['platform','mapping_accuracy','readout','new_assumptions','verdict']);w.writerows([['programmable mixed-signal categorical sampler + physical RNG','highest/designable arbitrary logits','direct 12-label','calibrated RNG+controller','SELECTED_DESIGN_ONLY'],['FB-MOSFET Potts p-bit array','promising probabilistic multi-state hardware; exact q12/A7 not demonstrated','one-hot state','extension/calibration to q12 arbitrary logits','CANDIDATE_COMPONENT'],['coupled CMOS ring oscillators','energy/optimization mapping stronger than exact Gibbs-rate mapping','phase state','sampling/rate proof missing','NO_GO_FOR_DYNAMIC_FIN_CURRENTLY']])
write('PHYS-006/NONCONCLUSIONS.md',nonc+'''\n- Nie stwierdzamy, że gotowe urządzenie q=12 spełniające kontrakt istnieje. To DESIGN_ONLY; niczego nie kupiono i nie kontaktowano laboratoriów.\n''')
write('PHYS-006/REPLAY.md','''# REPLAY\nBrak hardware replay. Matematyczny N=2 target: `python ../run_smallN_campaign.py`. Źródła literaturowe są jawnie podane w DEVICE_MAPPING.md.''')

# PHYS007
contr=pred['countermodel_metrics']; minTV=pr['minimum_primary_TV']; eps=0.003; robust=minTV-2*eps; ideal=max(contr[n]['M_equal_prior_bound_5pct'] for n in pr['primary_countermodels'])
radius=robust/2; a=radius/6; conservative=math.ceil(math.log(24/0.05)/(2*a*a))
r7={'execution_status':'COMPLETED_PREREGISTRATION_ONLY','scientific_verdict':'PASS_CONDITIONAL_STATIC_FALSIFICATION_DESIGN','claim_class':'DESIGN_ONLY+EXACT_MODEL_PREDICTIONS','validation_data_opened':False,'primary_min_TV':minTV,'per_model_TV_calibration_envelope_required':eps,'robust_primary_gap_after_two_envelopes':robust,'ideal_worst_primary_equal_prior_Chernoff_M_5pct':ideal,'very_conservative_12bin_union_bound_M_for_empirical_TV_radius':conservative,'money_cost':'OPEN until real platform throughput/cost is measured'}
dump('PHYS-007/RESULTS.json',r7)
# update prereg
pr.update({'per_model_total_TV_uncertainty_acceptance':eps,'robust_primary_gap':robust,'primary_test_rule':'pre-registered multinomial likelihood-ratio/model check only after calibration envelope passes','closest_models_below_robust_resolution':['PERT_36_5P','LEAK_K12_2P_TRACE']});dump('PHYS-007/PREREGISTRATION.json',pr)
dump('PHYS-007/RAW_RECORD_SCHEMA.json',{'fields':['run_id','dataset_role(calibration|design|validation)','timestamp_or_cycle','anchor_label','sampled_label','g_command','logit_vector_id','device_id','attempt_index','valid_one_hot','readout_quality_flag'],'validation_rule':'raw validation records immutable; no retuning g/A/readout after opening'})
write('PHYS-007/POWER_AND_ERROR.md',f'''# PHYS-007 — power and error budget\n\nPierwszy test zamrożono jako **N=2 conditional pair-difference histogram**, `g=3`, plus `g=0` negative control. To test statyczny: anchor `i=0`, próbka `j` ma dokładną predykcję `softmax[(g/2)A[0,j]]`. Nie używa bariery 327, protokołu 333 ani task335.\n\nPrimary finite alternatives: full Potts, flat P7, PERT_34_10P, PERT_35_10P. Minimalny idealny TV w tej klasie wynosi {minTV:.6f}. Przed walidacją każdy model+aparatura musi wejść w zamrożoną całkowitą kopertę kalibracyjną TV <= {eps:.3f}; po odjęciu dwóch kopert pozostaje dodatnia separacja >= {robust:.6f}. Dwie bliższe sondy (5% perturbacja k3/k6 i 2% leakage k1/k2) są tylko sensitivity checks, bo ich kontrast nie przechodzi tego budżetu.\n\nPrzy idealnie znanych rozkładach najgorszy primary Chernoff equal-prior bound <=5% wymaga około {ideal:,} niezależnych prób. To nie jest uniwersalna kontrola type-I. Skrajnie konserwatywny 12-bin union/Hoeffding bound dla empirycznego promienia TV odpowiadającego połowie robust gap wynosi ~{conservative:,} prób; pokazuje koszt pełnej distribution-free ostrożności. Realny plan ma użyć prerejestrowanego multinomial likelihood testu/parametric calibration po uzyskaniu rzeczywistej macierzy błędu odczytu.\n\nFinite-g approximation error dla N=2 pair law = 0 wewnątrz modelu. Niepewność `g` ma być kalibrowana niezależnie; cel ±0.2% przesuwa FIN prediction o maks. TV≈{pr['max_TV_shift_from_g_tolerance']:.6f} i musi mieścić się w całej kopercie 0.003.\n\nJeśli próbki pochodzą z jednej trajektorii, należy użyć ESS z autokorelacji; liczby powyżej zakładają niezależne reset/anchor cycles. Koszt pieniężny pozostaje OPEN, bo brak istniejącej zatwierdzonej platformy i throughput.\n''')
write('PHYS-007/REPORT.md',f'''# PHYS-007 REPORT\n\n**execution_status:** COMPLETED_PREREGISTRATION_ONLY  \n**scientific_verdict:** PASS_CONDITIONAL_STATIC_FALSIFICATION_DESIGN\n\nIstnieje pierwszy test, który nie zależy od metastability ani wyboru kinetyki: dokładny N=2 pair histogram. Primary countermodels są skończenie oddalone i pozostają rozdzielone po jawnej kopercie kalibracyjnej. Bliższe perturbacje nie są rozstrzygalne przy tej samej kopercie i zostały uczciwie zdegradowane do sensitivity-only zamiast retunowania testu. Żadne validation records nie zostały otwarte.\n''')
write('PHYS-007/NONCONCLUSIONS.md',nonc+'''\n- Sample counts są warunkowe na kalibrację i model błędu; nie są obietnicą wydajności nieistniejącego urządzenia.\n''')
write('PHYS-007/REPLAY.md','''# REPLAY\n`python ../run_smallN_campaign.py`\nNie ma danych validation. Nie wolno po ich przyszłym otwarciu zmieniać g=3, listy primary modeli ani koperty TV i nazywać tego tym samym testem.''')

# Inputs hashes per task
for d in [f'PHYS-{i:03d}' for i in range(1,8)]:
 p=ROOT/d/'INPUTS.sha256'
 refs=[ROOT/'inputs/FIN_PHYSICS_BRIDGE_MASTER_ROADMAP_20260929.md',ROOT/'inputs/HANDOFF_327_CURRENT.md',ROOT/'inputs/EPISTEMIC_LEDGER.md',ROOT/'run_smallN_campaign.py',ROOT/'PHYS-002/safe_seed.py']
 lines=[]
 for x in refs:
  lines.append(hashlib.sha256(x.read_bytes()).hexdigest()+'  '+str(x.relative_to(ROOT)))
 p.write_text('\n'.join(lines)+'\n')
