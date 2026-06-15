# P2777/S1727 symmetry-source selector/geometry audit

Status: `P2777_SYMMETRY_SOURCE_SELECTOR_GEOMETRY_AUDIT_NO_CLOSURE`

## Selector result
- orbit_of_directed_candidates=[1, 5, 7, 11]
- aut_invariant_singleton_selector_count=0
- symmetry_selects_orientation=False

## Geometry result
- all_pair_geometries_vertex_transitive=True
- max_automorphism_count_geometries=['torus_4x4']
- max_symmetry_selects_unique_geometry_on_pair=True

## Decision
Aut(Z12) symmetry cannot select +1 over -1.  Maximal automorphism count distinguishes the P2774 pair, but no strict law says to maximize automorphism count, no global strict graph class is audited, and no K/L_total variational coupling is exported.

## Recommendation
Do not claim selector or geometry closure from symmetry alone.  The next honest move is exactly one of: supply a strict source law that justifies a concrete symmetry functional such as maximal automorphism count and audit it over a declared 16-point/strict graph class, or supply a symmetry-breaking/chiral/pseudoscalar source that couples to the Aut(Z12) +1/-1 selector orbit.  Without one of those, preserve the P2697-P2777 no-selector/no-canonical-geometry/no-closure certificate.
