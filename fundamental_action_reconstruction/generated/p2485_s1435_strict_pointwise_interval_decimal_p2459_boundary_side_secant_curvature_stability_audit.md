# P2485/S1435 strict pointwise interval-Decimal P2459 boundary-side secant-curvature stability audit

Status: `PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_BOUNDARY_SIDE_SECANT_CURVATURE_STABILITY_AUDIT_NO_ROOT_WINDOW_NO_ANALYTIC_MONOTONICITY_CONVEXITY_NO_COVERAGE_INCREASE_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM`

## Boundary-side secant-curvature stability replay

Root-window-side boundary anchor: `4.058893978208217`.
Extended dyadic levels replayed: `64`.
Consecutive secant margins computed: `63`.
Consecutive secant-margin drifts computed: `62`.
Weakest extended Decimal separation: `7.54264363563016891325322230162706736783639817340995842904150075037988681078925444528958E-7`.
Tightest extended Decimal separation: `7.55236823851199790060667180175843591955048868586205588266515382435750345611443402981575E-7`.
Minimum consecutive positive lower-bound delta: `1.05434366776983891253621142342315916577851264569647597556699308452E-28`.
Minimum positive secant margin per removed width: `0.00248949820709680863101222591389247853875128857960824090250486742062874525155201631666176`.
Maximum positive secant margin per removed width: `0.00248949859905106794768353021303754319184229957320452770625378463619456441719577479479296`.
Secant margin absolute spread: `3.9195425931667130429914506465309101099359628680374891721556581916564375847813120E-10`.
Secant margin relative spread: `1.57443077564517987016880136538230675852319177826751343212552763958146459708102417341700033E-7`.
Minimum positive secant-margin drift: `8.499154928925678674706408368676348365327204979739680036618240E-29`.
Maximum positive secant-margin drift: `1.9597710960206805158181059396915797412734468281467716799612125220282626326702080E-10`.
Tightest extended level: `64`.
Tightest width fraction of P2482 leftmost subcell: `1/18446744073709551616`.
All extended rows exclude zero: `True`.
All secant margins positive: `True`.
All secant-margin drifts positive: `True`.

## Coverage budget

P2485 fresh Decimal evaluation rows (not a P2459 coverage count): `64`.
P2485 consecutive secant margins (not a P2459 coverage count): `63`.
P2485 consecutive margin drifts (not a P2459 coverage count): `62`.
Diagnostic row ratio against inherited P2459 universe (not a coverage fraction): `64/99846`.
New P2459 unreplayed cells added by P2485: `0`.
New P2459 unreplayed-cell scope against inherited P2459 universe: `0/99846`.
P2485 refines one inherited P2456 covered-boundary-chain cell: `True`.
Full complement replay exported by P2485: `False`.

## Plain-language progress note

This packet extends the one-sided boundary ladder to 64 dyadic levels and checks both positive finite secant margins and positive consecutive margin drifts. The finite gains behave coherently on the checked ladder, but this is still diagnostic arithmetic outside the excluded root window; it does not prove root-window, analytic monotonicity/convexity, selector/source, or ToE closure.

## Hard limits / negative controls

This is a finite boundary-side secant-curvature stability audit outside the excluded root window.  It is not a P2459 coverage increase, root-window exclusion theorem, full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity/convexity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.
