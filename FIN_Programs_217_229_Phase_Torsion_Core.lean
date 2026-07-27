/-!
FIN Programs 217--229: dependency-free phase-torsion finite core.

The analytic automatic-continuity and pi-transcendence arguments remain
paper mathematics.  This file machine-checks only the finite order-eight
endomorphism orbit and the reduced rational numerator/denominator datum of
the strict angle.
-/

def torsionImage8 (n : Nat) : Nat := n % 8

def imageOrbit8 : List Nat := (List.range 8).map torsionImage8

theorem image_orbit_is_complete :
    imageOrbit8 = [0, 1, 2, 3, 4, 5, 6, 7] := by
  decide

theorem every_image_is_torsion_slot :
    imageOrbit8.all (fun n => n < 8) = true := by
  decide

theorem strict_angle_fraction_is_reduced :
    Nat.gcd 743 4000 = 1 := by
  decide

#eval imageOrbit8
