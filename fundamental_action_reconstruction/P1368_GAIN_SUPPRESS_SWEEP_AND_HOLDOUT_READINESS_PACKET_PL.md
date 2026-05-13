# P1368 Gain/Suppress Sweep and Holdout-Readiness Packet (PL)

Status: `P1368_EXECUTED_SWEEP_AND_READINESS_NO_FALSE_PASS`
As of: `2026-05-12`
Artifacts:
- `generated/p1368_gain_suppress_sweep_summary.json`
- `generated/p1367b_postrun_governance_loop_summary.json`

## Cel

Wykonać kolejny krok po `P1367`:

1. sweep parametrów `gain/suppress`,
2. wybrać najlepszy punkt wg jawnego celu (mean abs-z + regularizacja tuningu),
3. ocenić `holdout_ready`,
4. uruchomić pętlę governance po runie.

## Wynik sweep

1. Przebadano siatkę 25 konfiguracji.
2. Najlepszy punkt nadal ma bardzo wysoki błąd (`best_max_abs_z >> 5`).
3. `holdout_ready = false`.

Wniosek: obecna klasa mapowania jest za słaba numerycznie na obecnym zestawie,
więc holdout byłby przedwczesny.

## Wynik pętli governance

Po runie ponowiono pętlę `P1364/P1365 -> P1362/P1363`.
Docelowo (po dołączeniu nowych artefaktów) pętla ma status `PASS`.

## Decyzja profesorska

Następny uczciwy krok: `P1369_MODEL_CLASS_REFINEMENT_BEFORE_HOLDOUT`

1. zmienić klasę mapowania kernel->stałe (nie tylko gain/suppress),
2. utrzymać małą liczbę parametrów i jawny koszt regularizacji,
3. po każdej iteracji obowiązkowo: residuale -> gate -> artifact check.

## Dla laika

Sprawdziliśmy wiele ustawień i najlepsze nadal jest za słabe.
To nie znaczy, że idea jest zła — znaczy, że potrzebujemy lepszego „silnika liczenia” zanim zrobimy finalny egzamin.
