# P2898/S1848 single-defect relation import-boundary audit

Status: `P2898_SINGLE_DEFECT_RELATION_IMPORT_BOUNDARY_NO_CLOSURE`

## Finite single-defect audit
- torsor size: `12`
- circulant backgrounds: `4096`
- directed edge placements: `144`
- single-defect candidates: `589824`
- translation quotient defect orbits: `49152`
- defect orbit size histogram: `{'12': 49152}`
- source-neutral defect placements: `0`
- candidates requiring imported pointer: `589824`

## Boundary
P2898 audits the nearest non-circulant repair after P2897: one labelled directed-edge defect on an arbitrary circulant relation.  There are 4096 circulant backgrounds and 144 edge placements, hence 589824 single-defect candidates, but for each background and edge difference the 12 placements form a free translation orbit.  The defect can point only because its labelled placement is imported; quotient-level data supplies 0 source-neutral defect placements.  Therefore one-edge non-circulant relation defects do not export a nonimported basepoint/polarity law or source the 9/5 carrier origin.

## Recommendation
Do not promote one-edge defects, labelled defect locations, non-circulant perturbations, edge differences, circulant backgrounds, relation profiles, scalar scores, or unpointed free-torsor clocks to strict phase/origin sourcehood.  A next proof-grade move must either supply an explicit strict law that sources the defect placement itself and couples it to the 9/5 variational density, or pivot to a genuinely different typed object outside torsor/basepoint/scalar-score/relation/defect/support/orbit/Fourier/inventory constructions; otherwise preserve no-new-live-frontier.
