# P1998 / S948 — Strict Cutkosky channelwise backend-weight calibration witness

P1998 to następny uczciwy krok po P1997: zachowuje jawny state-sum kanał-po-kanale,
ale zastępuje ręczne wagi kanałów (`gg`,`gh`,`hh`) wagami wyliczonymi z backendowych
współczynników krzywiznowych B1.

- strict lane only,
- brak transferu roszczeń legacy,
- status pozostaje `OPEN_OBSTRUCTION_WITH_TRACE`.

Wynik: `Delta_opt(s)` nadal jest niezerowe (kontrolowane), ale diagnostyka jest mniej arbitralna,
bo kluczowe proporcje kanałowe pochodzą z backendu.
