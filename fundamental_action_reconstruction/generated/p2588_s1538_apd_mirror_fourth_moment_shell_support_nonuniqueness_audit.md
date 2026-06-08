# P2588/S1538 APD mirror fourth-moment shell support nonuniqueness audit

Status: `P2588_APD_MIRROR_FOURTH_MOMENT_SHELL_SUPPORT_NONUNIQUENESS_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE`

## Result

- Frontier atom under attack: `strict_dynamical_source_for_A_P_D`.
- Supports audited: `3`.
- Cardinality: `8`.
- Mirror center: `5.5`.
- Internal second/fourth shells: `14.0` / `98.0`.
- Shell support nonunique: `True`.
- Strict APD dynamic source exported: `False`.

## Interpretation

P2588 strengthens P2587 by adding a fixed fourth-moment shell on eight-point mirror supports.  The supports share endpoints, cardinality, reflection symmetry, and unweighted central second/fourth moments; fixed-support full moments recover positive weights, but selected off-node APD dynamics still changes with the support.

## Recommended next honest step

Do not promote endpoints, mirror symmetry, cardinality, or fixed second/fourth support shells into an APD support source. P2588 keeps all of those constraints and still obtains distinct supports with conditionally unique weights and support-dependent APD dynamics. The next honest step is a strict nadsoliton-derived support/density law or a genuine internal APD dynamics theorem.

## Negative controls

No fourth-moment shell support source, kurtosis-shell source, variance-shell source, support-selection source, finite-support source, positive-measure source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.

## Fingerprint

`1bc5432d608f63009f19e4bad9dbcdeeffcbdae37d41aedee66a823b5632b546`
