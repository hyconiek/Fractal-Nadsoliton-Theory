# GENERAL-QUARTIC-CLOSURE-12 — arbitrary retained-probe candidate

Status: **THEOREM CANDIDATE / FINITE-SUM VERIFIED ON GENERIC MIXED CONTROLS**.

The proved pure-k5 theorem suggests a more general formula.  For any retained
real probe phi in the seven-dimensional space define
  h = P_H(phi^2),
  k = P_H(phi * (A7 phi)).
The candidate leading coefficient is

  C(phi) = -12 <h,h>_u + 2 g <h,k>_u,

with an exactly vanishing g^2 coefficient.

For an eigenmode A7 phi=lambda phi this reduces to
  C=(2 g lambda-12)||P_H(phi^2)||_u^2,
and for the pure k=5 cosine/sine it is exactly the already proved
  lambda5^2(g lambda5-6)/144.

Two nontrivial mixed controls were recomputed from the full finite-sum generator
series, not from this formula:
- a generic seven-coordinate mixture: residual <=4.8e-15;
- phi=k3c+0.37 k4c: residual <=1.7e-14.
The g^2 coefficient was zero in both computations.

This is strong evidence but not yet a tensor-level proof for every phi.  Do not
promote it beyond theorem-candidate status until the quartic tensor identity is
exported symbolically or independently certified on a spanning polynomial basis.
