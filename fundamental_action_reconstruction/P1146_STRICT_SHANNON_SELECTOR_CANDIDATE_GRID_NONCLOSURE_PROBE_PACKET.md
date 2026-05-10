# P1146 Strict Shannon Selector Candidate Grid Nonclosure Probe Packet

Status: `P1146_EXECUTED_STRICT_SHANNON_SELECTOR_CANDIDATE_GRID_NONCLOSURE_PROBE_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Perform a concrete follow-up to `P1145`: run a reproducible strict-only probe
on `S_cand(d)` and decide whether even a minimal uniqueness-proxy signal is
present.

## Method

Probe object:

```text
S_cand(d) = exp(-(4 ln 2)d) * cos(omega*d+phi)/(1+beta*d^eta)
```

with strict tuple:

```text
omega=0.18575, phi=0.16250, beta=1.0, eta=1.8
```

Finite-grid audit on `d in [0,24]`, step `0.1`.

Proxy criterion (explicitly non-theorem):

1. selector-like uniqueness hint requires non-oscillatory single-sign profile,
2. oscillation/sign flips imply candidate alone cannot support strict-core
   selector closure.

## Result

From generated summary packet:

1. negative values are present,
2. sign changes are present,
3. uniqueness proxy fails,
4. `QW-2191` nonclosure audit remains `BLOCKED`.

## Artifact

Machine summary:

- `generated/p1146_strict_shannon_selector_candidate_grid_nonclosure_probe_summary.json`

Reproducible probe script:

- `p1146_strict_shannon_selector_candidate_grid_nonclosure_probe.py`

## Strict-rigor interpretation

`P1146` strengthens anti-false-pass discipline:

1. strict-side Shannon weighting is admissible as a candidate construction,
2. but this concrete candidate does not yet provide a selector-uniqueness
   closure signal,
3. therefore no strict-core closure claim is admissible.

## Non-claims

`P1146` does not claim:

1. theorem-level impossibility of all Shannon-weighted candidates,
2. discharge of `QW-2191`,
3. ToE closure.
