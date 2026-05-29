# Scratch strict-alpha Hebbian energy-refined kWTA selector probe

Status: energy-refined Hebbian kWTA closes the finite d5 basin conditionally; not a strict source theorem.

- Supports scanned: `792` five-active-node states on `Z_12`.
- Raw kWTA branch histogram: `{'1': 624, '2': 120, '3': 24, '10': 12, '21': 12}`.
- Energy-refined branch histogram: `{'1': 720, '2': 60, '4': 12}`.
- All retained energy-refined branches reach d5: `792/792`.
- Closure layers: `{'0': 12, '1': 615, '2': 156, '3': 9}`; max layer `3`.
- Deterministic replay reaches d5: `792/792` with histogram `{'0': 12, '1': 480, '2': 264, '3': 36}`.
- Minimum non-d5 learned-energy ascent over retained branches: `4`.
- Target replay: `q^5=256/243`, eta `9/5`.
- Honest read: the previous kWTA tie obstruction is closed only after adding an explicit learned-energy refinement premise.
- No false pass: no d5-source theorem, no QW-2191 discharge, no ToE closure.
