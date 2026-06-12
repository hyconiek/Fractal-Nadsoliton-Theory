# P2646/S1596 frozen-kernel compression-signature preregistration certificate

Status: `P2646_FROZEN_KERNEL_COMPRESSION_SIGNATURE_PREREGISTRATION_CERTIFICATE_NO_EMPIRICAL_CONFIRMATION_NO_FALSE_PASS`

## Content-first anti-duplication audit

This audit greps blind/frozen empirical content, compression signatures, phase/amplitude-invariant observables, legacy-vs-strict discriminators, and nonclosure guardrails before adding the certificate.

- `blind_frozen_empirical_content`: 979 hits
- `compression_signature_content`: 128 hits
- `phase_amplitude_invariant_content`: 959 hits
- `legacy_strict_discriminator_content`: 4638 hits
- `guardrail_nonclosure_content`: 13836 hits

## Locked discriminator

Observable family: `phase/amplitude-invariant denominator tail ratios R(a,b)=A(b)/A(a) and log-tail slopes sigma(a,b)=-log R/log(b/a)`.
On blind holdout data, the measured denominator/envelope tail ratio must lie below the preregistered midpoint threshold and the measured log-tail slope must lie above the preregistered midpoint threshold for the audited pairs, with no kernel retuning.
Fail the strict compression-signature claim if ratios are legacy-like, if thresholds pass only after refitting beta/eta/phase per dataset, or if a matched exponential/spline/control baseline removes the frozen-kernel advantage.

| near | far | strict tail ratio | legacy tail ratio | midpoint ratio threshold | strict slope | legacy slope | midpoint slope threshold |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 7 | `0.058474514372` | `0.943925233645` | `0.501199874009` | `1.459041812800` | `0.029656208766` | `0.744349010783` |
| 1 | 12 | `0.022572225844` | `0.901785714286` | `0.462178970065` | `1.525624744802` | `0.041602510285` | `0.783613627543` |
| 2 | 8 | `0.103696464996` | `0.944444444444` | `0.524070454720` | `1.634780690630` | `0.041231080096` | `0.838005885363` |
| 3 | 9 | `0.154611007685` | `0.944954128440` | `0.549782568063` | `1.699273677849` | `0.051536738287` | `0.875405208068` |

## Verdict

P2646 converts the P2644 compression/locality-bias reading and P2645 role-matrix route into a falsifiable, phase/amplitude-invariant holdout protocol. The exact strict-vs-legacy denominator gaps are large on the audited pairs, so the next empirical check is no longer a vague analogy: it is a locked tail-ratio/log-slope discriminator. But this is a preregistration certificate only; without blind data it does not prove beta sourcehood, role transfer, L_total, or ToE closure.

Decision: `PREREGISTERED_FROZEN_COMPRESSION_SIGNATURE_INTERFACE_ONLY__NO_EMPIRICAL_CONFIRMATION_NO_SOURCE_EXPORT`.
Full kernel now? `False`.
ToE closure now? `False`.

## Professorial closure path

- Run the blind frozen-kernel compression-signature holdout exactly as locked here, including phase masks and control baselines.
- In parallel, continue the target-independent beta-source theorem; empirical success cannot replace the source theorem but can prioritize it.
- If the holdout fails or requires retuning, demote strict compression to a useful model-level successor and do not reopen role-bearing L_total.
- If the holdout passes, rerun the role-transfer matrix only for the modified/compressed successor semantics, not for unchanged inverse-hierarchy.
- Only after beta source, selector/source QW-2191, role-transfer matrix, and empirical controls pass can frozen-kernel L_total be reconsidered.

## Next honest step

Do not add more node/offset/stride or amplitude rescues. The next honest step is either an actual blind holdout execution of the P2646 tail-ratio/log-slope compression signature, or a target-independent beta-source theorem; the strongest route is to run the holdout while separately proving the beta source.

## Negative exports

- `empirical_confirmation_claimed`: `False`
- `data_fit_performed`: `False`
- `kernel_retuning_allowed`: `False`
- `beta_source_exported`: `False`
- `unchanged_inverse_hierarchy_restored`: `False`
- `legacy_role_transfer_revalidated`: `False`
- `bridge_completion_exported`: `False`
- `q_w_2191_discharged`: `False`
- `role_bearing_ltotal_reenabled`: `False`
- `toe_closure_claimed`: `False`
