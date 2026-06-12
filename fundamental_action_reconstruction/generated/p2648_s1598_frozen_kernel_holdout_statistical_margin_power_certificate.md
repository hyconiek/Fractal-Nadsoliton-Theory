# P2648/S1598 frozen-kernel holdout statistical margin/power certificate

Status: `P2648_STATISTICAL_MARGIN_RULE_READY_BUT_NO_BLIND_DATA_NO_FALSE_PASS`

## Content-first anti-duplication audit

This packet greps statistical margin/power, frozen holdout decision, legacy-vs-strict controls, and nonclosure guardrails before adding the certificate.

- `statistical_margin_power_content`: 261 hits
- `frozen_holdout_decision_content`: 172 hits
- `legacy_strict_control_content`: 23 hits
- `source_nonclosure_content`: 13857 hits

## Familywise decision rule

Familywise alpha: `0.05` over `8` inequalities.
Bonferroni one-sided z: `2.497705474412`.
Ratio rule: `measured_tail_ratio + z*ratio_standard_error < midpoint_ratio_threshold for every preregistered pair`.
Slope rule: `measured_log_tail_slope - z*slope_standard_error > midpoint_slope_threshold for every preregistered pair`.

| near | far | strict ratio margin | strict slope margin | max ratio SE | max slope SE | legacy ratio margin | legacy slope margin |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 7 | `0.442725359636` | `0.714692802017` | `0.177252828315` | `0.286139742791` | `-0.442725359636` | `-0.714692802017` |
| 1 | 12 | `0.439606744221` | `0.742011117259` | `0.176004236178` | `0.297077107313` | `-0.439606744221` | `-0.742011117259` |
| 2 | 8 | `0.420373989724` | `0.796774805267` | `0.168304067085` | `0.319002706055` | `-0.420373989724` | `-0.796774805267` |
| 3 | 9 | `0.395171560378` | `0.823868469781` | `0.158213834428` | `0.329850127736` | `-0.395171560378` | `-0.823868469781` |

## Verdict

P2648 fixes the main remaining weakness of the P2647 harness: a nominal threshold pass is not enough.  The frozen holdout must clear every locked ratio/slope inequality after a familywise one-sided uncertainty penalty.  The strict synthetic prediction has positive margin budget, while the legacy fixture is already on the wrong side of every audited inequality.  Still, this is a statistical decision rule and power/margin certificate only; no blind data have been tested.

Decision: `STATISTICAL_MARGIN_RULE_READY_BUT_NO_BLIND_DATA__NO_EMPIRICAL_CONFIRMATION_NO_SOURCE_EXPORT`.
Real blind holdout executed? `False`.
Full kernel now? `False`.
ToE closure now? `False`.

## Professorial closure path

- Require the blind payload to include standard errors or certified uncertainty bounds for each tail ratio and log-tail slope.
- Apply the P2648 familywise rule before reading any physical meaning into a P2647 holdout run.
- If any inequality passes only without the uncertainty penalty, record empirical failure rather than a weak pass.
- If the holdout passes with uncertainty and frozen controls, rerun only the modified/compressed successor row of the role-transfer matrix.
- Keep beta-source and selector/source proof obligations separate from empirical margin success.

## Next honest step

Acquire a blinded payload with uncertainty estimates tight enough to meet the P2648 margins, then run P2647/P2648 exactly once without retuning; in parallel continue the target-independent beta-source theorem.

## Negative exports

- `blind_holdout_data_loaded`: `False`
- `empirical_confirmation_claimed`: `False`
- `statistical_significance_claimed_on_real_data`: `False`
- `synthetic_margin_promoted_to_evidence`: `False`
- `kernel_retuning_allowed`: `False`
- `beta_source_exported`: `False`
- `role_transfer_revalidated`: `False`
- `q_w_2191_discharged`: `False`
- `role_bearing_ltotal_reenabled`: `False`
- `toe_closure_claimed`: `False`
