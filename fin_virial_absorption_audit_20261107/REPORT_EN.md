# Virial and feedback absorption audit

Date: 2026-11-07. Parent: b5e77b1dcbd99bcc0827a759c234718f164f93c1.

The requested upstream audit accepts VD's virial-to-local-energy argument
and LF/ML's weighted absorption at their stated dependency boundary.
No blocking defect was found in the audited steps. This is a same-author
analytical audit, not independent review of the complete FIN research chain.

The new derivation explicitly retains spatial fluxes, initial and terminal
virial terms, both components of phase forcing, the QP initial chart shift,
and the even-channel trace defect. The historical rational constants are
valid. The actual maximum of a''/a for a=<x>^-6 is 147/16; the earlier
bound by 12 is conservative and needs no correction.

The forcing contribution epsilon^4 T is O_S(epsilon^2) at T=S epsilon^-2.
The cumulative estimate is converted to a weighted estimate by integration
by parts with its endpoint retained, never by differentiating an inequality.
A negative control demonstrates why the latter shortcut would fail.

The weighted feedback coefficient can be made at most 1/2 by choosing
mu once from the source constants. Returning to ordinary time costs
exp(2mu S). This proves a fixed-S estimate, not an all-time bound or a
physical exponential damping rate. Substituting CS's higher-order trace
and TP2's separate cubic trace estimate uses no modal damping conclusion.

Validation: 49 checks passed, comprising 31 pinned input hashes and 18
exact symbolic/constant/scaling controls. These check two virial flux
identities, the coercivity density transformation, rational margins,
weight curvature, temporal endpoints and critical-time exponents. They
are not a numerical nonlinear PDE evolution or formal theorem proof.

Next audit LE/QH's even-channel propagation and modulated causal tails,
including threshold regularity and switch-on transients used by RW.
SP, CM, CD and the exact spectral resolution remain inherited inputs;
independent review and empirical FIN validation remain open.
