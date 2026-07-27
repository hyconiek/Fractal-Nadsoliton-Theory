/-!
FIN Programs 178--190
Finite Dirichlet core, deliberately independent of the AFIS Boolean model.

This file uses only Lean core arithmetic.  It records a concrete four-cycle
certificate and the constant-vector null statement.  It does not formalize
matrix exponentials, complete positivity, or physical adequacy.

Release status: source authored; no machine-compilation claim is made unless
the release manifest explicitly reports a successful Lean toolchain run.
-/

def sq (x : Int) : Int := x * x

/- Twice the Dirichlet energy of the unit-row-sum four-cycle.
   The four undirected edges are (0,1), (1,2), (2,3), (3,0). -/
def twiceDirichlet4 (x0 x1 x2 x3 : Int) : Int :=
  sq (x0 - x1) + sq (x1 - x2) + sq (x2 - x3) + sq (x3 - x0)

theorem constant_vector_is_null (c : Int) :
    twiceDirichlet4 c c c c = 0 := by
  simp [twiceDirichlet4, sq]

theorem archived_test_vector :
    twiceDirichlet4 1 2 (-1) 3 = 30 := by
  decide

theorem archived_energy_nonnegative :
    0 ≤ twiceDirichlet4 1 2 (-1) 3 := by
  decide

/- The general positivity theorem is elementary over Int, but is intentionally
   left outside the release claim because a local Lean toolchain with the
   desired ordered-ring tactic library is not available. -/
