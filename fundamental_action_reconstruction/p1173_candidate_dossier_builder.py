#!/usr/bin/env python3
"""P1173 builds a formula-level markdown dossier for top candidate."""
from __future__ import annotations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
GEN=ROOT/'generated'


def main():
    top=json.loads((GEN/'p1170_neighbor_candidate_e2e_screen_summary.json').read_text())['rows'][0]
    cand=json.loads(Path(top['candidate']).read_text())
    safe=json.loads((GEN/'p1172_safe_operating_region_summary.json').read_text())

    omega=cand['omega_hint']; phi=cand['phi_hint']; beta=cand['beta_hint']; eta=cand['eta_hint']
    sigma=cand['sigma_hint']; kappa=cand['kappa_hint']

    md=f'''# P1173 Top Candidate Formula Dossier

Status: `P1173_EXECUTED_TOP_CANDIDATE_FORMULA_DOSSIER_NO_FALSE_PASS`
As of: `2026-05-10`

## Candidate identity

- Premise ID: `{cand['premise_id']}`
- Source file: `generated/{Path(top['candidate']).name}`

## Core formulas

Baseline strict kernel with asymmetric selector term:

```text
N_sigma(d) = cos(omega*d + phi) + sigma*(1 - exp(-kappa*d))
K_cand(d)  = exp(-alpha*d) * N_sigma(d) / (1 + beta*d^eta)
alpha      = 4 ln 2
```

Numerical parameter point:

```text
omega = {omega}
phi   = {phi}
beta  = {beta}
eta   = {eta}
sigma = {sigma}
kappa = {kappa}
```

## Operational safe region (from P1172)

```text
omega in [{safe['safe_bounds']['omega']['min']}, {safe['safe_bounds']['omega']['max']}]
phi   in [{safe['safe_bounds']['phi']['min']}, {safe['safe_bounds']['phi']['max']}]
sigma in [{safe['safe_bounds']['sigma']['min']}, {safe['safe_bounds']['sigma']['max']}]
kappa in [{safe['safe_bounds']['kappa']['min']}, {safe['safe_bounds']['kappa']['max']}]
```

## Strict-rigor interpretation

1. Candidate passes E2E workflow checks on current repo state.
2. Candidate has proxy-stable region rather than only a single isolated point.
3. This remains a `CANDIDATE_PASS_ONLY` methodological status.
4. No strict-core closure claim and no `QW-2191` discharge claim are made.
'''

    out=ROOT/'P1173_TOP_CANDIDATE_FORMULA_DOSSIER.md'
    out.write_text(md,encoding='utf-8')
    print(f"[P1173] wrote {out}")

if __name__=='__main__':
    main()
