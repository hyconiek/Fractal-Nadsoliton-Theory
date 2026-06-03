# P2482/S1432 strict pointwise interval-Decimal P2459 boundary-band weakest-cell dyadic-refinement replay audit

Status: `PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_BOUNDARY_BAND_WEAKEST_CELL_DYADIC_REFINEMENT_REPLAY_AUDIT_NO_ROOT_WINDOW_NO_COVERAGE_INCREASE_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM`

## Boundary-band weakest-cell dyadic refinement

Parent collar side: `p2456_right_boundary_band`, boundary-band index `0`.
Dyadic subcell replay rows: `128`.
Parent Decimal separation inherited from P2481: `5.06293652999729190722297749614506090018367516699953176293930422739711844797226919218097E-7`.
Minimum subcell Decimal separation: `7.53291903478976770323585651880162653212203235374303123554573410172935718772595565693193E-7`.
Maximum subcell Decimal separation: `7.90714975634730658078393789868540074477741278909008348705710749130107094132189697931211E-7`.
Minimum consecutive endpoint gap: `0E-15`.
Maximum consecutive endpoint gap: `0E-15`.
Minimum consecutive positive delta: `2.94653671210199220691274322451099717950475429995774284152710930827784065050168077444E-10`.
Minimum is leftmost subcell adjacent to root window: `True`.
Refined minimum lower bound exceeds parent interval lower bound: `True`.
Refined-minus-parent lower-bound delta: `2.46998250479247579601287902265656563193835718674349947260642987433223873975368646475096E-7`.
All subcells exclude zero: `True`.

## Coverage budget

P2482 fresh Decimal evaluation rows (not a P2459 coverage count): `128`.
Diagnostic row ratio against inherited P2459 universe (not a coverage fraction): `128/99846`.
New P2459 unreplayed cells added by P2482: `0`.
New P2459 unreplayed-cell scope against inherited P2459 universe: `0/99846`.
P2482 refines one inherited P2456 covered-boundary-chain cell: `True`.
Full complement replay exported by P2482: `False`.

## Plain-language progress note

This packet opens the weakest P2481 collar row, which is the first P2456 right-boundary-band cell next to the excluded root window. Its 128 dyadic diagnostic rows all exclude zero and increase left-to-right, but the weakest row is still the leftmost subcell on the root-window side. The refined lower bound is stronger than the coarse parent-cell lower bound. This is finite covered-boundary-cell refinement, not a proof inside the root window, not a P2459 coverage increase, and not a continuum or full-complement proof.

## Hard limits / negative controls

This is a finite dyadic refinement of one inherited P2456 covered-boundary-chain cell only.  It is not a P2459 coverage increase, root-window exclusion theorem, full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.
