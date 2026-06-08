# P2587/S1537 APD mirror second-moment shell support nonuniqueness audit

Status: `P2587_APD_MIRROR_SECOND_MOMENT_SHELL_SUPPORT_NONUNIQUENESS_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE`

## Result

- Frontier atom under attack: `strict_dynamical_source_for_A_P_D`.
- Supports audited: `3`.
- Cardinality: `6`.
- Mirror center: `5.5`.
- Internal second-moment shell: `25.0`.
- Shell support nonunique: `True`.
- Strict APD dynamic source exported: `False`.

## Interpretation

P2587 strengthens P2586 by adding a fixed second-moment shell on six-point mirror supports.  The supports share endpoints, cardinality, reflection symmetry, and unweighted central second moment; fixed-support full moments recover positive weights, but selected off-node APD dynamics still changes with the support.

## Recommended next honest step

Do not promote endpoints, mirror symmetry, cardinality, or a fixed second-moment support shell into an APD support source. P2587 keeps all of those constraints and still obtains distinct supports with conditionally unique weights and support-dependent APD dynamics. The next honest step is a strict nadsoliton-derived support/density law or a genuine internal APD dynamics theorem.

## Negative controls

No second-moment shell support source, variance-shell source, support-selection source, finite-support source, positive-measure source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.

## Fingerprint

`4a7ce2779910709d67e41669a4aa8a910a92686230328c25363400c870d18766`
