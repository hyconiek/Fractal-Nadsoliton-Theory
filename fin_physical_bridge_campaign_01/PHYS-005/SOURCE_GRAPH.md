# PHYS-005 — A7 source and necessity ledger

`supplied strict kernel W(d)` → `circulant Laplacian Fourier spectrum` → **selected** support `{3,4,5,6}` → feature matrix `X7` → `A7=X7 X7^T`.

Pierwsza i pogrubiona strzałka nie są obecnie niezależnie fizycznie wyprowadzone. D12, PSD i rank 7 ograniczają klasę, ale po ustaleniu aktywnego supportu pozostają cztery dodatnie wagi, a po trace-normalization trzy niezależne proporcje. Zatem symmetry/rank/trace nie wybiera konkretnego tuple FIN.

Gaussian-parent jest poprawnym mostem konstrukcyjnym:
`H_parent=1/2 x^T K x - N^-1/2 x^T sum_a c_{sigma_a} - h n0`.
Eliminacja Gaussowskiego x przez completion of square daje pair matrix proporcjonalną do `c_i^T K^-1 c_j`. Jednak dowolna PSD A ma taką faktoryzację; wybór K,c po to, by odtworzyć A7, jest engineered realization. Staje się źródłem predykcyjnym dopiero wtedy, gdy K i c są ustalone z niezależnej fizyki/pomiarów przed oglądaniem A7 fingerprint.

**Verdict: NO_NEW_SOURCE.** Zachowujemy najwyżej dwa falsyfikowalne atomy: niezależny mediator measurement i niezależną selection rule. Bez nich nie dodajemy nowej „zasady fundamentalnej”.
