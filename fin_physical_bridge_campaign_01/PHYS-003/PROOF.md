# PHYS-003 — static spectral fingerprint

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
