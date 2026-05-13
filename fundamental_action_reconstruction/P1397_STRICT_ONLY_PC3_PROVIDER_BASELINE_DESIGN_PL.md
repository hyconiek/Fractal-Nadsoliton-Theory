# P1397 — STRICT-ONLY PC3 PROVIDER BASELINE DESIGN (PL)

## Cel kroku
Wykonać pierwszy krok po `P1396`: zaprojektować nową klasę provider `PC3` w torze strict-only, tak by nie powtarzać cyklicznego strojenia `PC2` i utrzymać ruch w kierunku
`F_Nadsoliton ⇒ L_SM + L_GR` bez bridge do legacy.

## Repo grep audit (czy to już było robione?)
Przeprowadzono grep-audit pod kątem duplikacji kroku:
- brak wcześniejszego pakietu `P1397_STRICT_ONLY_PC3_PROVIDER_BASELINE_DESIGN`,
- brak wcześniejszego checkpointu `p1397_*`,
- `PC3` występuje tylko jako trigger/plan po `P1396`, nie jako gotowy baseline-run.

Wniosek: ten krok nie duplikuje istniejącego exportu wykonawczego.

## Zakres strict-only
- `legacy_bridge_used = false`
- brak transferu semantyki legacy -> strict
- kontynuacja po formalnym obstruction `PC2-DRIFT-v1`

## Baseline PC3 (v1)
- `provider_class_id = PC3_strict_selector_anchor_v1`
- `anchor_family = A3_noncyclic_selector_stabilizers`
- `inherits_from_pc2_loop = false`
- `epsilon_sign_v1 = 0.05` (zamrożone)
- `epsilon_drift_v1 = 0.04` (zamrożone)
- `first_run_contract = dual_metric_required`
  - metryka 1: `sign_flip_rate <= epsilon_sign_v1`
  - metryka 2: `selector_drift <= epsilon_drift_v1`

## Status kroku
- `PC3_BASELINE_STATUS := FROZEN_READY_FOR_FIRST_RUN`
- `L_B1_03_EXPORT_STATUS := NOT_EXPORTED`
- `B1_STATUS := OPEN`

## Decyzja profesorska (następny uczciwy krok)
Uruchomić `P1398_STRICT_ONLY_PC3_FIRST_DUAL_METRIC_RUN`: pierwszy run na zamrożonych progach i pełny werdykt PASS/FAIL bez retuningu kryteriów po fakcie.

## Omówienie dla laika
To jak start nowej generacji prototypu: poprzedni model miał oficjalnie udokumentowane ograniczenie, więc teraz budujemy nowy model od czystej bazy i testujemy go na tych samych, uczciwych progach jakości.
