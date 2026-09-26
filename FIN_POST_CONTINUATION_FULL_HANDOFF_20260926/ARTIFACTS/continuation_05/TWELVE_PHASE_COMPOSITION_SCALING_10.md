# GEOMETRY-TRANSITION-10 — 12-phase composition-controlled M=16 -> M=32 campaign

Status: **EXHAUSTIVE FINITE-SIZE M=16 COMPOSITION CERTIFICATE + EXHAUSTIVE M=32 CURVATURE AT ALPHA=0.425 + SELECTED M=32 PEAK CERTIFICATES; NO THERMODYNAMIC THEOREM**.

The 495 four-doubled M=16 compositions form 29 D12 orbits.  `fin12_m16_peaks.cpp`
uses exact automatic differentiation through the count-vector partition recursion
and brackets the heat-capacity-like maximum for every orbit.

Weighted over all 495 compositions, M=16 has
- beta_peak mean ~0.63503, standard deviation ~0.03192;
- full orbit range beta_peak = 0.5955796 ... 0.7664867;
- C_H,peak weighted mean ~7.0749, range 6.5025 ... 7.9377.

Thus the 12-phase finite-M pseudocritical point is strongly composition dependent.

Replicate each composition once to M=32 (counts 4 on the same four labels and 2
on the other eight).  `fin12_m32_jet.cpp` differentiates log Z directly; no
finite-difference derivative is used.  All 29 D12 orbits, hence all 495
compositions, were evaluated at beta=0.2125, corresponding to alpha=2 beta=0.425.

At this common scaled point:
- weighted mean C_H(M=16, beta=0.425) = 5.26889;
- weighted mean C_H(M=32, beta=0.2125) = 50.48048;
- M=32 full composition range = 48.84887 ... 52.32899;
- the composition-by-composition correlation between C_16 and C_32 is ~0.9623;
- the C_32/C_16 ratio has weighted mean ~9.590, range ~8.906 ... 10.158.

This robust curvature does **not** imply alpha-collapse of the peak.  Exact 3rd
and 4th derivatives show, for example:
- mask 15: M=16 beta*=0.6971175, while M=32 beta*=0.2136003 (alpha*=0.4272005);
- mask 195: M=16 beta*=0.7664867, M=32 local Newton beta*~0.2134194;
- mask 585: M=16 beta*=0.5955796, M=32 local Newton beta*~0.2156113;
- mask 275: M=16 beta*=0.5961539, M=32 local Newton beta*~0.2151137.

Therefore the simple 8-phase scaling variable alpha=beta*n is **not** a valid
collapse law for the 12-phase replicated 2/1 -> 4/2 composition sequence.
Interestingly, selected M=32 peak locations are much more tightly clustered
than their M=16 ancestors; this is finite-size evidence that composition
sensitivity may contract under replication, not a limit theorem.

For mask 15 the energy cumulants at the peaks change from
- M=16: skewness ~0.7297, Binder-like U4 ~ -0.1049;
- M=32: skewness ~0.2854, Binder-like U4 ~ +0.1966.

The root-composition order histogram also does not support a clean two-state
mixture.  At the M=16 peak mask 15 has effective macrobin count ~1.87 and ~77%
weight in the maximally segregated bin.  At the M=32 peak it has effective
macrobin count ~7.24 and three local maxima (probabilities ~0.198, 0.260, 0.158).
Mask 585 gives the same qualitative M=32 pattern with effective macrobin count
~7.49 and maxima ~0.188, 0.238, 0.171.  The strongest local barrier between the
q2=20 and q2=24 sectors is only order log-probability 1.2--1.7.

Conclusion: 12-phase M=32 has a large, composition-robust heat-capacity feature,
but its relation to M=16 is not the 8-phase alpha-collapse and its root landscape
is multimodal rather than a clean binary phase mixture.  M=64 or an analytic
fixed-composition limit is still required before any thermodynamic claim.
