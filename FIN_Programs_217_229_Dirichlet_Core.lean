/-!
FIN Programs 217--229: dependency-free finite Dirichlet certificate.

This source deliberately uses only Lean core.  It machine-checks two
independent integer formulas for the C4 Dirichlet form on archived witnesses,
the constant null mode, and nonnegativity of the archived energy.  It is not
the general Mathlib theorem.
-/

def sq217 (x : Int) : Int := x * x

def edgeDirichlet4 (x0 x1 x2 x3 : Int) : Int :=
  sq217 (x0 - x1) + sq217 (x1 - x2) +
  sq217 (x2 - x3) + sq217 (x3 - x0)

def laplacianPairing4 (x0 x1 x2 x3 : Int) : Int :=
  2 * (sq217 x0 + sq217 x1 + sq217 x2 + sq217 x3) -
  2 * (x0*x1 + x1*x2 + x2*x3 + x3*x0)

theorem constant_edge_null (c : Int) :
    edgeDirichlet4 c c c c = 0 := by
  simp [edgeDirichlet4, sq217]

theorem constant_laplacian_null (c : Int) :
    laplacianPairing4 c c c c = 0 := by
  simp [laplacianPairing4, sq217]

theorem archived_dirichlet_identity :
    laplacianPairing4 1 2 (-1) 3 = edgeDirichlet4 1 2 (-1) 3 := by
  decide

theorem archived_energy_value :
    edgeDirichlet4 1 2 (-1) 3 = 30 := by
  decide

theorem archived_energy_nonnegative :
    0 ≤ edgeDirichlet4 1 2 (-1) 3 := by
  decide

#eval edgeDirichlet4 1 2 (-1) 3
