# P2778/S1728 max-symmetry 16-node geometry source audit

Status: `P2778_MAX_SYMMETRY_16NODE_GEOMETRY_SOURCE_AUDIT_NO_CLOSURE`

## Finite class result
- candidate_count=19
- max_automorphism_count=4096
- max_geometry_names=['circulant_pm1_pm7', 'circulant_pm3_pm5']
- torus_4x4_automorphism_count=384
- max_symmetry_selects_torus_4x4=False

## Decision
The supplied max-automorphism candidate law fails on the declared 16-node class: it does not select torus_4x4 and is not unique on labeled circulant candidates.  No K/L_total variational coupling is exported.

## Recommendation
Do not use maximal automorphism count as the geometry source.  The next honest move is either to test a different sourced symmetry functional with a declared target and quotient rule, or pivot from symmetry to a genuine strict metric/variational source such as a spectral action with fixed target spectrum and K/L_total coupling.  Otherwise preserve the P2697-P2778 no-canonical-geometry/no-closure certificate.
