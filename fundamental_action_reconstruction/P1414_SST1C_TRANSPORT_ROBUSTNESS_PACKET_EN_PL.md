# P1414 SST-1C Transport Robustness Packet (EN/PL)

Status: `P1414_EXECUTED_STRICT_TRANSPORT_TEST_NO_FALSE_PASS`
As of: `2026-05-13`

## Professorial decision

After `SST-1B` replay convergence, the next honest step is transport robustness,
not promotional narrative.

## Execution

- Script: `p1414_sst1c_transport_robustness_checkpoint.py`
- Summary artifact: `generated/p1414_sst1c_transport_robustness_summary.json`
- Obstruction export (if fail): `generated/p1414_selector_obstruction_v1.json`

## Result

- `verdict = FAIL_STRICT_TRANSPORT_DRIFT`
- strict promotion toward `L_SM + L_GR` remains blocked.

This is an honest scientific result: reproducibility improved (`SST-1B`),
but robustness still fails under transport perturbations.

## Lay explanation (PL)

Dla laika: dwa zespoły policzyły podobnie (to był dobry znak),
ale kiedy lekko zmieniamy warunki testu, wynik „ucieka” za mocno.
To znaczy: model nie jest jeszcze stabilny i nie wolno ogłaszać sukcesu.

## Recommendation

Freeze this failure as formal obstruction and move to a **new noncyclic provider class**
for selector stabilization (strict-only, no legacy bridge).
