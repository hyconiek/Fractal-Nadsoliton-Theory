# P2484/S1434 strict pointwise interval-Decimal P2459 boundary-side dyadic secant margin replay audit

Status: `PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_BOUNDARY_SIDE_DYADIC_SECANT_MARGIN_REPLAY_AUDIT_NO_ROOT_WINDOW_NO_ANALYTIC_MONOTONICITY_NO_COVERAGE_INCREASE_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM`

## Boundary-side dyadic secant margin replay

Root-window-side boundary anchor: `4.058893978208217`.
Secant levels replayed: `32`.
Consecutive secant margins computed: `31`.
Weakest secant Decimal separation: `7.54264363563016891325322230162706736783639817340995842904150075037988681078925444528958E-7`.
Tightest secant Decimal separation: `7.55236823850746952903590999826300668543652162351399169169616736413082116024936917615770E-7`.
Minimum consecutive positive lower-bound delta: `4.52837157181614672053475381162428464555836141700267701276176963351545449120E-19`.
Minimum positive secant margin per removed width: `0.00248949820709680863101222591389247853875128857960824090250486742062874525155201631666176`.
Maximum positive secant margin per removed width: `0.0024894985990510675826476056642992188405847185850079636626256352275244521664358809337856`.
Tightest secant level: `32`.
Tightest width fraction of P2482 leftmost subcell: `1/4294967296`.
Tightest-minus-P2482 lower-bound delta: `1.944920371770182580005347946138015331448926977096045615043326240146397252341351922577E-9`.
All secant rows exclude zero: `True`.
All consecutive secant margins positive: `True`.

## Coverage budget

P2484 fresh Decimal evaluation rows (not a P2459 coverage count): `32`.
P2484 consecutive secant margins (not a P2459 coverage count): `31`.
Diagnostic row ratio against inherited P2459 universe (not a coverage fraction): `32/99846`.
New P2459 unreplayed cells added by P2484: `0`.
New P2459 unreplayed-cell scope against inherited P2459 universe: `0/99846`.
P2484 refines one inherited P2456 covered-boundary-chain cell: `True`.
Full complement replay exported by P2484: `False`.

## Plain-language progress note

This packet replays 32 one-sided dyadic contractions next to the excluded root window and computes 31 finite secant margins from the lower-bound gains. Every checked row excludes zero and every normalized finite gain is positive. The result is useful boundary-side margin evidence, but it is still finite diagnostic arithmetic outside the root window and does not prove root-window, analytic, selector/source, or ToE closure.

## Hard limits / negative controls

This is a finite boundary-side dyadic secant-margin audit outside the excluded root window.  It is not a P2459 coverage increase, root-window exclusion theorem, full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.
