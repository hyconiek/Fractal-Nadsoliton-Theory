# P2488/S1438 strict pointwise interval-Decimal P2459 boundary-handoff collar monotonicity lemma certificate

Status: `PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_BOUNDARY_HANDOFF_COLLAR_MONOTONICITY_LEMMA_CERTIFICATE_PROOF_COMPRESSION_NO_NEW_REPLAY_NO_ROOT_WINDOW_NO_GLOBAL_ANALYTIC_MONOTONICITY_NO_COVERAGE_INCREASE_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM`

## Finite collar monotonicity lemma

Total collar rows compressed: `134`.
Boundary-band rows: `6`.
Dyadic subcell rows: `128`.
All lemma preconditions met: `True`.
Finite piecewise monotone-increasing collar lemma exported: `True`.
Finite positive collar zero-exclusion lemma exported: `True`.
Minimum row amplitude lower bound: `5.06293652999729190722297749614506090018367516699953176293930422739711844797226919218097E-7`.
Minimum derivative lower bound: `0.000376525332868308529985112453876382902868708201825128622178302466106607764170990234336364200`.
Total checked collar width: `0.0006999999999980000000`.
Derivative transport lower gain over checked collar: `2.63670284865778124097316884096847344624475419545110405715616914967156550658511793601823408E-7`.

## Coverage budget

New Decimal replay rows in P2488: `0`.
Reused P2487 derivative rows (not a P2459 coverage count): `134`.
Reused diagnostic row ratio against inherited P2459 universe (not a coverage fraction): `134/99846`.
New P2459 unreplayed cells added by P2488: `0`.
New P2459 unreplayed-cell scope against inherited P2459 universe: `0/99846`.
P2488 is proof compression, not new replay: `True`.

## Plain-language progress note

This packet compresses the P2487 derivative sweep into a finite collar lemma. The checked collar rows are exactly adjacent, positive, and have positive derivative intervals, so the checked collar is piecewise monotone increasing and zero-excluding. This is a collar-local proof-compression step and does not prove the excluded root window, global analytic monotonicity, selector/source closure, or ToE closure.

## Hard limits / negative controls

This is a proof-compression lemma for the finite checked boundary-handoff collar outside the excluded root window.  It is not a P2459 coverage increase, root-window exclusion theorem, full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, global analytic monotonicity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.
