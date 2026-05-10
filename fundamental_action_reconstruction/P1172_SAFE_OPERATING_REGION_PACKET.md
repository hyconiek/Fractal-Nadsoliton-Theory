# P1172 Safe Operating Region Packet

Status: `P1172_EXECUTED_SAFE_OPERATING_REGION_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Wykonać rekomendowany krok po `P1171`: wyznaczyć praktyczny bezpieczny region
parametrów dla top kandydata, gdzie proxy pozostaje stabilne.

## Professor-level decision

Buduję 5x5x5x5 siatkę perturbacji (`625` punktów) wokół top kandydata i
wyodrębniam punkty spełniające:

- `negative_count=0`,
- `sign_change_count=0`

na domenie `[0,72]`.

## Result

- `stable_points = 408/625` (`stable_fraction = 0.6528`),
- wyeksportowano `safe_bounds` dla `(omega,phi,sigma,kappa)`.

To daje pierwszą operacyjną mapę regionu roboczego zamiast pojedynczego punktu.

## Artifacts

- script:
  `p1172_safe_operating_region.py`
- summary:
  `generated/p1172_safe_operating_region_summary.json`

## Honest boundary

`P1172` jest mapą regionu proxy-stabilnego, nie dowodem closure ani discharge
`QW-2191`.
