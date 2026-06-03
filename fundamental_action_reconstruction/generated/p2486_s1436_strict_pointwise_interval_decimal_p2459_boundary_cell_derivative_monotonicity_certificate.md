# P2486/S1436 strict pointwise interval-Decimal P2459 boundary-cell derivative monotonicity certificate

Status: `PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_BOUNDARY_CELL_DERIVATIVE_MONOTONICITY_CERTIFICATE_LOCAL_ONE_CELL_NO_ROOT_WINDOW_NO_GLOBAL_ANALYTIC_MONOTONICITY_NO_COVERAGE_INCREASE_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM`

## One-cell derivative monotonicity certificate

Cell: `[4.058893978208217, 4.058894759458217]`.
Cell width: `7.81250000E-7`.
Amplitude interval: `{'lo': '7.53291903478976770323585651880162653212203235374303123554573410172935718772595565693193E-7', 'hi': '7.57476429222637067872945729008398201002106890052745618431309569592786331842692301709842E-7'}`.
Amplitude interval excludes zero: `True`.
Amplitude separation from zero: `7.53291903478976770323585651880162653212203235374303123554573410172935718772595565693193E-7`.
Derivative interval: `{'lo': '0.000377193506523256731912425579069501927841700655278864552965861600570292173765314779757572356', 'hi': '0.000377200389077054648538724623302131833331034154137127973526400599938835968580526332909181426'}`.
Derivative interval positive on entire cell: `True`.
Derivative separation from zero: `0.000377193506523256731912425579069501927841700655278864552965861600570292173765314779757572356`.
Derivative relative width: `0.0000182467451822166154869028154406431084083235363031612385194216564706294726568638872579318811`.
Local finite-interval monotone-increasing witness exported: `True`.

## Coverage budget

P2486 derivative interval rows (not a P2459 coverage count): `1`.
Diagnostic row ratio against inherited P2459 universe (not a coverage fraction): `1/99846`.
New P2459 unreplayed cells added by P2486: `0`.
New P2459 unreplayed-cell scope against inherited P2459 universe: `0/99846`.
P2486 refines one inherited P2456 covered-boundary-chain cell: `True`.
Full complement replay exported by P2486: `False`.

## Plain-language progress note

This packet checks the derivative on the whole already-inherited boundary-side subcell instead of only replaying smaller dyadic samples. The projection amplitude is positive on the cell, and its interval derivative is also strictly positive, so the cell has a local monotone-increasing witness. The certificate remains local to one covered boundary-chain cell outside the excluded root window and does not prove global analytic, selector/source, or ToE closure.

## Hard limits / negative controls

This is a local one-cell derivative monotonicity certificate outside the excluded root window.  It is not a P2459 coverage increase, root-window exclusion theorem, full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, global analytic monotonicity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.
