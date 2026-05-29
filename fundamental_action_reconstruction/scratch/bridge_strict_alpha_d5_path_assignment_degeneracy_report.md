# Scratch strict-alpha d5 path assignment degeneracy probe

Status: assignment degeneracy certificate after distance-5 support selection; no assignment selector theorem.

- Reference ordered support `[0, 7, 2, 9, 4]` has all adjacent cyclic distances equal to `5`.
- One-support assignment table: `10` balanced assignments with transition histogram `{'1': 2, '3': 4, '2': 3, '4': 1}`.
- Orbit scan: `24` ordered supports and `240` assignment rows.
- Path-smooth minimizers: energy `1`, count `48`, unique high-node sets `12`.
- Target replay: `q^5=256/243`, eta residual `0.000e+00`.
- Honest read: path smoothing narrows assignments to endpoint blocks, but does not choose endpoint, translate, orientation, or source phase.
- No false pass: no strict assignment theorem, no QW-2191 discharge, no ToE closure.
