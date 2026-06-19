# P2907/S1857 joint origin-sign source theorem candidate gate

Status: `P2907_JOINT_ORIGIN_SIGN_SOURCE_THEOREM_CANDIDATE_GATE_NO_STRICT_EXPORT`

## Finite gate
- Xi family before joint postulate: `24`
- Xi family after joint postulate: `1`
- selected origin/sign: `(0, 1)`
- selected defect edge: `[0, 5]`
- nonzero local derivatives: `1`
- strict provenance exported: `False`
- unit-bearing L_total coupling exported: `False`

## Boundary
P2907 constructs the missing joint origin-and-sign object as an explicit theorem-postulate J_{0,+}.  The finite gate confirms that it collapses the 24 Xi alternatives to Xi_{0,+}, selects D=(0,5), and gives one symbolic local derivative U_9_5.  This is proof-readiness only: J_{0,+} is still imported, U_9_5 is still symbolic, and no strict nadsoliton provenance or unit-bearing nonproxy L_total coupling is exported.

## Recommendation
Audit/prove strict provenance for the joint J_{0,+} theorem itself: a nadsoliton-derived rule must compute both origin 0 and sign + without importing them.  If that cannot be supplied, do not add more postulated J variants; pivot outside Xi/defect-placement or preserve no-new-live-frontier.  Only after provenance may the next lane lift U_9_5 to a unit-bearing L_total coupling theorem.
