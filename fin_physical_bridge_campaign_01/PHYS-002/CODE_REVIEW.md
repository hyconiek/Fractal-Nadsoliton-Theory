# PHYS-002 — reproduction and seed-code audit

## Verdict
`PASS_SCOPED_RECONSTRUCTED_CORE__ORIGINAL_FIN_SON_UNAVAILABLE`. Repozytoryjne, śledzone źródła wystarczają do niezależnego odtworzenia małego rdzenia, ale dwa historyczne pliki `FIN son` nie są dostępne w bazowym commicie, więc nie twierdzę, że zostały naprawione bajt-w-bajt.

## Fail-closed checks
A7: rank=7, trace=15.2613460079, centering residual=2.005e-15, diagonal spread=2.220e-16, min eigenvalue=-3.876e-16. Wszystkie trzy N=2 generatory przechodzą row-sum/stationarity/detailed-balance bez symetryzowania błędnego wejścia.

Klasyfikacja modów używa projektorów/orbit C12, a nie etykiety pierwszego wektora z degenerowanego multipletu. Przy N=2,g=0 wszystkie 12 sektorów są obecne; maksymalna eigenpair residual jest 3.026e-15.

Dla N=3, heat-bath, `g=G_FROZEN` sektor k=4 daje `rho=0.13143978619564556`, różnica od fixture 2.220e-16.

`G_FROZEN` zgadza się z R09 `g_bal` do ~2e-15 i ma pochodzenie operacyjne (balans kanałów bariery), ale nie jest stałą fizyczną. `V_d4` nie jest używane jako dowód globalny. Aktualny scoped intake nadal wymaga dokładnych enclosure wejść, LP i scalar minimum przed promocją `Gamma=B4`.
