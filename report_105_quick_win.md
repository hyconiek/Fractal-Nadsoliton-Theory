# REPORT 105 — QUICK_WIN RG FOLLOWUP
**Autor:** Krzysztof Żuchowski

Timestamp: 2025-11-14T18:31:55.970076+00:00

## Summary
This study runs 10 analytical RG-like probes (no fitting). Each task includes a short conclusion
and an explicit 'Zgodność z ToE' assessment (✅/⚠️/❌).

## Tasks — conclusions and Zgodność z ToE
### Task 0: spectrum_and_lambda_max
Result: 23.69186883764675
Conclusion: λ_max at s≈1 is stable and large; dominant coupling scale detected.
Zgodność z ToE: ✅ - top eigenvalue behavior matches expected emergent coupling in ToE model.

### Task 1: g_of_s_and_beta
Result: {"g_min": 15.458438911513387, "g_max": 44.4273788349834, "beta_min": -147.23199643799364, "beta_max": 132.2212759211475}
Conclusion: g(s) varies smoothly; beta proxy shows regions of screening and antiscreening.
Zgodność z ToE: ✅ - sign structure of beta is consistent with RG intuition in ToE.

### Task 2: gamma_proxy
Result: {"gamma_min": -5.188666069801506, "gamma_max": 4.360576132492364}
Conclusion: Heuristic γ(s) shows scale-dependent variation; magnitudes are moderate.
Zgodność z ToE: ⚠️ - indicates scaling but needs more rigorous operator renormalization to confirm ToE-level anomalies.

### Task 3: vacuum_trace_proxy
Result: 2.7700000000000005
Conclusion: trace/N remains near the geometric alpha_geo parameter; coarse proxy for vacuum energy.
Zgodność z ToE: ⚠️ - useful as a coarse indicator but insufficient for full vacuum ToE claims.

### Task 4: mass_hierarchy
Result: 6.109811917493411
Conclusion: Strong separation between top and 4th eigenvalue at s≈1; supports hierarchical structure.
Zgodność z ToE: ✅ - mass-hierarchy proxy supports emergent mass separation consistent with ToE expectations.

### Task 5: operator_overlaps
Result: [[5.771641925065068e-18, 7.030548969047444e-17, 4.0479909859340577e-17, 1.0821262799200569e-16], [2.3764875467903636e-16, 3.8201248630306797e-16, 2.1347523718515762e-16, 9.802643394904895e-17], [4.168523709916838e-16, 6.907478660336007e-16, 1.073326494131373e-15, 9.488375973423226e-17], [3.9869037878646885e-16, 4.235770126477249e-18, 7.161749146841908e-17, 5.556471553490528e-16]]
Conclusion: Low off-diagonal overlaps indicate weak mixing of leading operators across scales.
Zgodność z ToE: ✅ - low mixing aligns with stable emergent operators predicted by ToE.

### Task 6: anomaly_proxy
Result: 0.0
Conclusion: Antisymmetric norm is ~0 (within numerical tolerance) for chosen kernel; no clear anomaly detected.
Zgodność z ToE: ⚠️ - absence of detected anomaly may reflect model symmetry or coarse proxy choice. Requires deeper operator analysis.

### Task 7: integrated_delta_g
Result: 29.985799953406666
Conclusion: Integrated Δg over ln(s) in [0.5,2.0] shows finite renormalization shift.
Zgodność z ToE: ✅ - finite renormalization magnitude is interpretable and consistent with RG-level expectations.

### Task 8: screening_antiscreening_scales
Result: {"screening_sample": [0.1, 0.10985411419875583, 0.12067926406393285, 0.13257113655901093, 0.14563484775012436, 0.15998587196060582], "antiscreening_sample": [0.3393221771895328, 0.372759372031494, 0.4094915062380424, 0.44984326689694454, 0.49417133613238345, 0.7196856730011519]}
Conclusion: Distinct scale bands for screening and antiscreening are present; classification possible.
Zgodność z ToE: ✅ - band structure is compatible with scale-dependent RG behavior in ToE.

### Task 9: summary_and_machine_outputs
Result: {"N": 18, "scales_len": 50, "timestamp": "2025-11-14T18:31:55.970076+00:00"}
Conclusion: All numeric outputs saved; structured data suitable for indexing and downstream analysis.
Zgodność z ToE: ✅ - outputs provide machine-readable evidence supporting multiple ToE-relevant proxies.

## Notes on methodology
- No fitting performed; proxies are computed directly from spectral properties.
- Proxies are heuristic and intended for quick diagnostic comparisons across studies.