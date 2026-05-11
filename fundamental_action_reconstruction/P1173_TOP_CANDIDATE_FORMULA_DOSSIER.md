# P1173 Top Candidate Formula Dossier

Status: `P1173_EXECUTED_TOP_CANDIDATE_FORMULA_DOSSIER_NO_FALSE_PASS`
As of: `2026-05-10`

## Candidate identity

- Premise ID: `P1170_NEIGHBOR_1`
- Source file: `generated/p1170_neighbor_1.json`

## Core formulas

Baseline strict kernel with asymmetric selector term:

```text
N_sigma(d) = cos(omega*d + phi) + sigma*(1 - exp(-kappa*d))
K_cand(d)  = exp(-alpha*d) * N_sigma(d) / (1 + beta*d^eta)
alpha      = 4 ln 2
```

Numerical parameter point:

```text
omega = 0.15
phi   = 0.09
beta  = 0.95
eta   = 1.75
sigma = 1.2
kappa = 0.12
```

## Operational safe region (from P1172)

```text
omega in [0.13, 0.16999999999999998]
phi   in [0.06, 0.12]
sigma in [1.0999999999999999, 1.4]
kappa in [0.07999999999999999, 0.16]
```

## Strict-rigor interpretation

1. Candidate passes E2E workflow checks on current repo state.
2. Candidate has proxy-stable region rather than only a single isolated point.
3. This remains a `CANDIDATE_PASS_ONLY` methodological status.
4. No strict-core closure claim and no `QW-2191` discharge claim are made.
