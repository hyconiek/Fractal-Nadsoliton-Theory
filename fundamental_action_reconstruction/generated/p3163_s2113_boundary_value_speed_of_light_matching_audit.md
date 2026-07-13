# P3163/S2113 boundary-value speed-of-light matching audit

Status: `P3163_BOUNDARY_VALUE_SPEED_OF_LIGHT_MATCHING_UNDERDETERMINATION_AUDIT`

## Theorem steps
- `C1` `receiver_definition`: Let v_hat be any nonzero dimensionless boundary velocity receiver from the strict model.
- `C2` `dimension_lift`: A physical velocity has dimension U_length/U_time, so v_phys = v_hat*(U_length/U_time).
- `C3` `fit_degeneracy`: For any target c>0, choosing U_length/U_time=c/v_hat gives v_phys=c.
- `C4` `no_closure`: Therefore numerical agreement with c is not a source theorem unless U_length and U_time are independently strict-sourced and a Lorentz/light embedding is exported.

## Finite certificate
- `boundary_value_candidates`: `10`
- `speed_of_light_fit_rows`: `10`
- `numerically_fit_nonzero_candidates`: `9`
- `scale_degeneracy_rows`: `225`
- `candidate_gate_rows`: `80`
- `accepted_physical_c_sources`: `0`
- `c_SI_m_per_s`: `299792458.0`

## Finite theorem
`P3163_T1_boundary_c_fit_is_scale_underdetermined`: Boundary/limiting model values can be used as receiver diagnostics for physical constants, but any nonzero dimensionless v_hat can match c by choosing U_length/U_time=c/v_hat.  Since current artifacts do not export independent U_length, U_time, Lorentz/spacetime embedding, or observed-light source, no boundary-value fit to c is a strict dimension source.

## Decision
P3163 constructs the requested boundary-value-to-c matching audit.  Nine nonzero boundary receivers can numerically fit c after choosing a unit ratio, and the zero tail cannot; all successful fits are scale choices rather than strict source laws.

## Recommendation
If this boundary-value route is pursued, the next proof-grade object must be a two-axis unit theorem U_LT: independently source U_length and U_time (or a Lorentzian metric/light-cone embedding) before comparing to c.  Otherwise preserve P3162's S_+ unit-source frontier or pivot to Lambda_origin_source_localizer for phase-origin work.
