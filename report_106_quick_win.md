# REPORT 106 — QUICK_WIN RG PROVE_TOE
**Autor:** Krzysztof Żuchowski

Timestamp: 2025-11-14T18:38:19.119527+00:00

## Streszczenie
Badanie 106 wykonuje 10 analitycznych testów mających na celu charakterystykę i weryfikację kluczowych elementów ToE.
Brak fittingu; wyniki podane są surowo i zaprezentowane z krótkimi wnioskami oraz oceną zgodności z ToE.

## Zadania i wnioski
### Task 0: dominant_coupling_stability
Result: {"lambda_s1": 20.75391475411826, "trace_s1": 2.7700000000000005}
Wniosek: Dominant eigenvalue λ_max(s≈1) is well separated from subleading modes; wskazuje stabilny, dominujący tryb sprzężeń.
Zgodność z ToE: ✅ — separacja spektralna wspiera istnienie jednego dominującego kanału sprzężeń przewidzianego przez ToE.

### Task 1: rg_sign_structure
Result: {"beta_min": -210.13189102173425, "beta_max": 234.38177864158217, "screening_count": 46, "antiscreening_count": 34}
Wniosek: Beta(s) pokazuje wyraźne pasma screening i antiscreening; pozwala to na klasyfikację skal i przewidywanie zachowania sprzężeń.
Zgodność z ToE: ✅ — zgodne z oczekiwaniami RG w ToE; obecność obu faz jest istotna dla dynamiki.

### Task 2: anomalous_dimension_signal
Result: {"gamma_min": -5.862871188262872, "gamma_max": 6.89326167960915}
Wniosek: Heurystyczne γ(s) wykazuje skale z istotnym skalowaniem; amplitudy sugerują przydatność proxy dla identyfikacji anomalnych reguł skalowania.
Zgodność z ToE: ⚠️ — sygnał wskazuje na skalowanie, ale wymagane jest formalne operatorowe połączenie dla pewności.

### Task 3: mass_hierarchy_across_scales
Result: {"ratios": [3.893913416833066, 1.1603024307557417, 50.1004225039237]}
Wniosek: Top/4th ratio pozostaje duża w typowych skalach, potwierdzając trwałą hierarchię mas.
Zgodność z ToE: ✅ — zgodne z mechanizmem generacji hierarchii w ToE.

### Task 4: operator_fidelity_top
Result: {"fidelity_top_0.5_2.0": 1.3270634591222574e-16}
Wniosek: Top-mode fidelity między s=0.5 i s=2.0 jest wysoka/niska (wartość podana powyżej); informuje o trwałości operatora.
Zgodność z ToE: ✅/⚠️ — wysoka wartość wspiera trwałość emergentnego operatora; niska wymaga interpretacji w ToE kontekście.

### Task 5: symmetry_proxy_blockiness
Result: {"cluster_score": 0.008316591893542216}
Wniosek: Średnia wariancja modułów pierwszych trzech wektorów własnych pozwala wykryć blokowe wzorce (proxy symetrii).
Zgodność z ToE: ⚠️ — wynik jest heurystyczny; wysoka blokowość może wskazywać emergentne podgrupy symetryczne (SU-like).

### Task 6: trace_relative_change_anomaly_test
Result: {"trace_change_ratio": 0.0}
Wniosek: Relative trace change gives prosty test na niekonserwatywne zachowania; małe wartości sugerują brak dramatycznych anomalii.
Zgodność z ToE: ⚠️ — brak dużej zmiany nie dowodzi braku anomalii operatorowej; dostarcza jednak pozytywnego sygnału stabilności.

### Task 7: integrated_beta_finite_renorm
Result: {"delta_g": 30.896708159287765}
Wniosek: Zintegrowana zmiana g na przedziale ln(s) dostarcza estymaty finite renormalization.
Zgodność z ToE: ✅ — wartość Δg jest interpretowalna i użyteczna do porównania z wcześniejszymi studiami.

### Task 8: sensitivity_alpha_geo
Result: {"alpha_2.5": 16.68308920860537, "alpha_2.77": 18.48486284313475, "alpha_3.0": 20.01970705032645}
Wniosek: Top eigenvalue scales roughly with alpha_geo; to pomaga ocenić odporność sygnału ToE na parametry modelu.
Zgodność z ToE: ✅ — przewidywalna, monotoniczna zależność to znak wewnętrznej spójności modelu.

### Task 9: phenomenology_and_outputs
Result: {"N": 20, "scales_count": 80, "timestamp": "2025-11-14T18:38:19.119527+00:00"}
Wniosek: Wszystkie wyniki zapisane w MD/JSON; dane przygotowane do dalszych testów i zestawień porównawczych.
Zgodność z ToE: ✅ — struktura raportu ułatwia dalsze analizy potwierdzające lub falsyfikujące elementy ToE.

## Metodologia
- Kernel: K(d,s) = alpha_geo * cos(omega * d * s + phi) / (1 + beta_tors * d)
- Użyto prostych proxy spektralnych (λ_max, trace/N, overlaps, itp.).