# P1807 S757 Strict Unified Runpack To State-Vector Lock Sync Packet (PL)

Status: `P1807_EXECUTED_STRICT_UNIFIED_RUNPACK_TO_STATE_VECTOR_LOCK_SYNC_PACKET_NO_FALSE_PASS`

## Cel

Po `P1806` domknąć brakujący krok operacyjny:
spójna synchronizacja `run-pack verdict -> state-vector theorem locks`
bez ręcznej interpretacji i bez fałszywych promocji gate status.

## Zakres

Wejście:

- `generated/p1806_s756_strict_unified_nonproxy_ea_eh_elg_residual_runpack_contract_checkpoint.json`

Wyjście:

- zaktualizowany artefakt state-vector lock sync dla `TG1/TG2/TG3`.

## Reguły synchronizacji (twarde)

1. Jeśli `TG1_BW != PASS_ZERO`:
   - `SV6_BW_global = OPEN_LOCKED_BY_UNIFIED_NONPROXY_RESIDUAL`,
   - `SV7_BRST_global = OPEN_LOCKED_BY_TG1`,
   - `SV8_Cutkosky_global = OPEN_LOCKED_BY_TG2`.

2. Jeśli `TG1_BW = PASS_ZERO` ale brak BRST witness:
   - `SV6_BW_global = PASS_ZERO`,
   - `SV7_BRST_global = OPEN_LOCKED_BY_NILPOTENCY_WITNESS`,
   - `SV8_Cutkosky_global = OPEN_LOCKED_BY_TG2`.

3. Zakaz ręcznego ustawienia `SV7/SV8 = PASS` bez odpowiadających witness-packów (`P1801/P1802`).

## Co zostało dowiedzione

1. Przejście z kontraktu wykonawczego (`P1806`) do state-vector jest teraz deterministyczne.
2. Lock chain `BW -> BRST -> CUT` jest utrzymany jako mechanizm automatyczny.

## Co pozostaje OPEN

1. Faktyczne `PASS_ZERO` dla trzech residuali w `P1806`.
2. BRST nilpotency witness (`P1801`), Cutkosky witness (`P1802`).

## Ryzyka false-pass

1. Ręczna edycja state-vector z pominięciem checkpointu `P1806`.
2. Semantyczne mylenie `UNLOCKED` z `PASS`.

## Następny uczciwy krok

Uruchomić realny nonproxy run-pack na wspólnym freeze,
a następnie przepuścić wynik przez ten sync, bez manual override.
