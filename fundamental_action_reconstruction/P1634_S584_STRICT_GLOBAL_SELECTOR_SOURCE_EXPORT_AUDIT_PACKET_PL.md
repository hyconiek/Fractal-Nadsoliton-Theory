# P1634 / S584 — Strict global selector source export audit (from existing atlas objects)

## Cel
Nie powielać badań: wykorzystać istniejące obiekty atlas/cocycle w repo i sprawdzić,
czy dają pełny eksport `E_selector_internal_source_full_domain` wymagany przez strict-core closure.

## Tor
`K_strict -> coeff -> L_total -> EOM -> local invertibility -> overlap compatibility -> selector source export audit`.

## Zasada
Brak false-pass: jeśli obiekty są lane-scoped lub bez globalnego domknięcia operator-level,
status pozostaje `KEEP_OPEN`.
