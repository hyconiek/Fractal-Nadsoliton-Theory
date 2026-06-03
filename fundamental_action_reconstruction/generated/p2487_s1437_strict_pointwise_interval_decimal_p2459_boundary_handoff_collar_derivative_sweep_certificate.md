# P2487/S1437 strict pointwise interval-Decimal P2459 boundary-handoff collar derivative sweep certificate

Status: `PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_BOUNDARY_HANDOFF_COLLAR_DERIVATIVE_SWEEP_CERTIFICATE_FINITE_PIECEWISE_MONOTONE_NO_ROOT_WINDOW_NO_GLOBAL_ANALYTIC_MONOTONICITY_NO_COVERAGE_INCREASE_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM`

## Boundary-handoff collar derivative sweep

Collar interval: `[4.058893978208217, 4.059593978208215]`.
Boundary-band derivative rows: `6`.
Dyadic subcell derivative rows: `128`.
Total derivative rows: `134`.
All amplitude intervals positive: `True`.
All derivative intervals positive: `True`.
All rows have local monotone-increasing witness: `True`.
All consecutive rows exactly adjacent: `True`.
All consecutive amplitude separations strictly increase: `True`.
Minimum amplitude interval separation: `5.06293652999729190722297749614506090018367516699953176293930422739711844797226919218097E-7`.
Maximum amplitude interval separation: `0.000001016933134107930802946894248688882194434615869381378311750866075781669598796178650793367`.
Minimum derivative lower bound: `0.000376525332868308529985112453876382902868708201825128622178302466106607764170990234336364200`.
Maximum derivative upper bound: `0.000377616529328645225995460343173670720689331574189524448573106162816365865352822669835717100`.
Maximum derivative relative width: `0.00233865230401300794842798299251727475741392520928573678723712733796446445768650377313620574`.

## Coverage budget

P2487 derivative interval rows (not a P2459 coverage count): `134`.
Diagnostic row ratio against inherited P2459 universe (not a coverage fraction): `134/99846`.
New P2459 unreplayed cells added by P2487: `0`.
New P2459 unreplayed-cell scope against inherited P2459 universe: `0/99846`.
P2487 reuses P2481 collar rows without new P2459 coverage: `True`.
Full complement replay exported by P2487: `False`.

## Plain-language progress note

This packet checks the derivative on every diagnostic row of the P2481 handoff collar. All 134 rows have positive projection amplitude and positive projection-amplitude derivative, and the rows remain exactly adjacent. This provides finite piecewise monotonicity evidence for the checked seam, but it is still outside the excluded root window and does not prove global analytic, selector/source, or ToE closure.

## Hard limits / negative controls

This is a finite piecewise derivative monotonicity certificate for the checked boundary-handoff collar outside the excluded root window.  It is not a P2459 coverage increase, root-window exclusion theorem, full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, global analytic monotonicity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.
