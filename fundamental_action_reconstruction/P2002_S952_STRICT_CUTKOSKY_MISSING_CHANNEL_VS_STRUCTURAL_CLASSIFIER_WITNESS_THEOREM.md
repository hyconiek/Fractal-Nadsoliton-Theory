# P2002 / S952 — Strict Cutkosky missing-channel vs structural classifier witness

P2002 wykonuje następny uczciwy krok po P2001: porównuje normy `Delta_opt`
pomiędzy P2000 i P2001, aby skategoryzować trend po dodaniu kanału `gx`.

Klasyfikator:
- spadek norm -> `MISSING_CHANNEL_PRESSURE_SUPPORTED`,
- wzrost norm -> `STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED`,
- brak wyraźnego trendu -> `MIXED_OR_INCONCLUSIVE`.

To jest diagnostyka porównawcza strict-lane, a nie theorem-grade closure.
