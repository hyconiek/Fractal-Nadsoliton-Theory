/-!
Finite Lean 4 certificate for the logical independence witnesses used in
FIN Program 167.  This file uses only Lean core: no Mathlib theorem is needed.

The Booleans record whether an AFIS model supplies axioms A0,...,A5.  The
certificate does not formalize the analytic Dirichlet-form theorem; it
formalizes the finite capability/independence layer and keeps that scope
explicit.
-/

structure AFISFlags where
  a0 : Bool
  a1 : Bool
  a2 : Bool
  a3 : Bool
  a4 : Bool
  a5 : Bool
  deriving Repr, DecidableEq

def full : AFISFlags := ⟨true, true, true, true, true, true⟩

def withoutA0 : AFISFlags := { full with a0 := false }
def withoutA1 : AFISFlags := { full with a1 := false }
def withoutA2 : AFISFlags := { full with a2 := false }
def withoutA3 : AFISFlags := { full with a3 := false }
def withoutA4 : AFISFlags := { full with a4 := false }
def withoutA5 : AFISFlags := { full with a5 := false }

def dualDynamics (m : AFISFlags) : Bool :=
  m.a0

def targetState (m : AFISFlags) : Bool :=
  m.a0 && m.a1

def genericOutcome (m : AFISFlags) : Bool :=
  m.a0 && m.a1 && m.a3

def signedOutcome (m : AFISFlags) : Bool :=
  m.a0 && m.a1 && m.a2 && m.a3

def calibratedSignedTest (m : AFISFlags) : Bool :=
  m.a0 && m.a1 && m.a2 && m.a3 && m.a4

def typedCompletion (m : AFISFlags) : Bool :=
  m.a0 && m.a1 && m.a5

def roleAuditEligible (m : AFISFlags) : Bool :=
  m.a0 && m.a1 && m.a2 && m.a3 && m.a4 && m.a5

theorem a0_independent :
    withoutA0.a0 ≠ full.a0 ∧
    withoutA0.a1 = full.a1 ∧ withoutA0.a2 = full.a2 ∧
    withoutA0.a3 = full.a3 ∧ withoutA0.a4 = full.a4 ∧
    withoutA0.a5 = full.a5 := by decide

theorem a1_independent :
    withoutA1.a1 ≠ full.a1 ∧
    withoutA1.a0 = full.a0 ∧ withoutA1.a2 = full.a2 ∧
    withoutA1.a3 = full.a3 ∧ withoutA1.a4 = full.a4 ∧
    withoutA1.a5 = full.a5 := by decide

theorem a2_independent :
    withoutA2.a2 ≠ full.a2 ∧
    withoutA2.a0 = full.a0 ∧ withoutA2.a1 = full.a1 ∧
    withoutA2.a3 = full.a3 ∧ withoutA2.a4 = full.a4 ∧
    withoutA2.a5 = full.a5 := by decide

theorem a3_independent :
    withoutA3.a3 ≠ full.a3 ∧
    withoutA3.a0 = full.a0 ∧ withoutA3.a1 = full.a1 ∧
    withoutA3.a2 = full.a2 ∧ withoutA3.a4 = full.a4 ∧
    withoutA3.a5 = full.a5 := by decide

theorem a4_independent :
    withoutA4.a4 ≠ full.a4 ∧
    withoutA4.a0 = full.a0 ∧ withoutA4.a1 = full.a1 ∧
    withoutA4.a2 = full.a2 ∧ withoutA4.a3 = full.a3 ∧
    withoutA4.a5 = full.a5 := by decide

theorem a5_independent :
    withoutA5.a5 ≠ full.a5 ∧
    withoutA5.a0 = full.a0 ∧ withoutA5.a1 = full.a1 ∧
    withoutA5.a2 = full.a2 ∧ withoutA5.a3 = full.a3 ∧
    withoutA5.a4 = full.a4 := by decide

theorem full_has_every_capability :
    dualDynamics full = true ∧
    targetState full = true ∧
    genericOutcome full = true ∧
    signedOutcome full = true ∧
    calibratedSignedTest full = true ∧
    typedCompletion full = true ∧
    roleAuditEligible full = true := by decide

theorem removing_a2_loses_signed_not_generic :
    genericOutcome withoutA2 = true ∧
    signedOutcome withoutA2 = false := by decide

theorem removing_a4_loses_calibration_not_signed_outcome :
    signedOutcome withoutA4 = true ∧
    calibratedSignedTest withoutA4 = false := by decide

theorem removing_a5_loses_completion_not_calibrated_test :
    calibratedSignedTest withoutA5 = true ∧
    typedCompletion withoutA5 = false ∧
    roleAuditEligible withoutA5 = false := by decide

#eval full
#eval calibratedSignedTest full
#eval roleAuditEligible withoutA5
