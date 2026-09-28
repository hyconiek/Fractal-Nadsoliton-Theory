# PHYS-004 — countermodels and kinetic robustness

**execution_status:** COMPLETED  
**scientific_verdict:** PASS_STATIC_SPECIFICITY__DYNAMIC_NOT_A7_ONLY

Przed porównaniem zamrożono trace A i użyto tego samego g. Zbudowano full centered Potts, flat rank-7 P7, trzy perturbacje trace-preserving oraz osobny 2% leakage do k=1,2.

Statyczny N=2 pair histogram przy g=3 odróżnia FIN m.in. od full Potts (TV≈0.096999), flat P7 (TV≈0.011787) i dwóch 10% perturbacji aktywnych wag (TV≈0.0128–0.0148). Najbliższe zaplanowane sondy 5%/2% mają TV≈0.00442/0.00563 i są znacznie trudniejsze.

Dynamika nie jest A7-only. Po jawnej normalizacji zegara `rate_k4(g=0)=1`, przy tym samym A7 i g=G_FROZEN sektor k=4 ma: heat-bath=0.265087, Metropolis=0.056810, Barker=0.102037. Każda kinetyka zachowuje tę samą pi, ale relaksuje inaczej.

Wniosek: statyka może testować kształt A7 niezależnie od kinetyki; R_k wymaga osobnego, jawnego kontraktu generatora.
