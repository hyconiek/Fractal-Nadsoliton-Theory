# P1416 PC7 First Transport Run Packet (EN/PL)

Status: `P1416_EXECUTED_PC7_RUN1_NO_FALSE_PASS`
As of: `2026-05-13`

## Professorial decision

Run the first PC7 transport test with strict threshold-only adjudication.

## Execution

- Script: `p1416_pc7_first_transport_run_checkpoint.py`
- Summary: `generated/p1416_pc7_first_transport_run_summary.json`
- Obstruction (if fail): `generated/p1416_pc7_margin_obstruction_v1.json`

## Result

- Boundary drift and replay gap passed,
- but selector margin failed preregistered minimum.
- Final verdict: `FAIL_PC7_MARGIN`.

Scientific interpretation:

PC7 improves stability but does not yet secure selector separation strength.
So promotion toward `L_SM + L_GR` remains blocked.

## Lay explanation (PL)

Po ludzku: nowa wersja PC7 jest stabilniejsza niż wcześniej,
ale „różnica jakości” między najlepszym wyborem a alternatywami jest nadal za mała.
To jak sygnał, który jest już czystszy, ale wciąż za słaby, żeby zaufać mu bez wątpliwości.

## Recommendation

Pivot to `PC8` noncyclic selector-margin amplifier and rerun transport + replay under the same strict thresholds.
