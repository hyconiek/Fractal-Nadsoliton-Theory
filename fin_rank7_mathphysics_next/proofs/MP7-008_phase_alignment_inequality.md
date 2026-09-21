# MP7-008 — exact phase-alignment inequality

Use the coefficients of MP7-007.  Since `c_{-m}=c_m`, the series is real and

`Z(0;b)-Z(phi;b) = sum_m c_m [1-cos(m.phi)] >= 0`

for every finite `a3,a4,a5>=0` and every `b>=0`.

Hence

`Z(phi;b) <= Z(0;b)`.

Equality holds precisely when `m.phi in 2*pi*Z` for every frequency `m` with strictly positive coefficient.  This equality condition is classified in MP7-009/010.

For `b<0`, label translation is exact.  If

`T_l: phi_k -> phi_k+2*pi*k*l/12,  b -> (-1)^l b`,

then relabeling `j -> j+l` gives `Z(T_l(phi,b))=Z(phi,b)`.  Choosing odd `l` converts a negative alternating amplitude to positive before applying the inequality.  No signed `sinh(b)` expansion is used directly.

The finite-sum diagnostic `results/MP7-007_010_phase_alignment_diagnostics.json` checks 236 adversarial/random phase points and the negative-`b` translation identity; it is a regression check, not the proof.
