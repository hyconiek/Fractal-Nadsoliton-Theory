# P1999 / S949 — Strict Cutkosky backend-kappa calibrated channel solver witness

P1999 wykonuje następny uczciwy krok po P1998: kalibruje z backendu nie tylko
wagi kanałów, ale także parametry odpowiedzi kanałowej (`kappa_i`).

- strict lane only,
- brak transferu roszczeń legacy,
- status pozostaje `OPEN_OBSTRUCTION_WITH_TRACE`.

Wynik: `Delta_opt(s)` pozostaje niezerowe, ale model proxy ma mniejszą arbitralność,
bo oba poziomy kalibracji kanałowej (wagi + kappa) są backend-zasilane.
