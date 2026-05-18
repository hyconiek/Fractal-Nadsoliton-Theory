# P2017 S967 Strict Cutkosky Quadrature-Ansatz Provenance Witness Theorem

Status: `P2017_CORRECTED_QUADRATURE_ANSATZ_PROVENANCE_GATE_NO_FALSE_PASS`
As of: `2026-05-18`

## Goal

Correct and sharpen the P2017 strict-lane Cutkosky/unitarity step.  P2017 is not
allowed to claim that a full backend loop-amplitude theorem has been exported.
It exports a useful finite quadrature **candidate** and, crucially, a provenance
gate explaining why the candidate remains below theorem-grade backend closure.

The channel remains

```text
graviton -> gauge_gauge.
```

## Strict kernel used

The only kernel used is the operational strict working kernel:

```text
K_strict(d) = cos(omega*d + phi)/(1 + beta*d^eta)
```

with the strict tuple

```text
omega = 0.18575, phi = 0.16250, beta = 1.0, eta = 1.8.
```

No legacy physical-role transfer is used.

## Candidate construction

For each exported `P2004` energy sample `s` and channel

```text
i ∈ {gg, gh, hh, gx},
```

P2017 defines a finite quadrature chart with `x ∈ [0,1]`, channel distance

```text
d_i(s,x) = sqrt(s) * (1 + x + offset_i),
offset_gg = 0, offset_gh = 1/4, offset_hh = 1/2, offset_gx = 3/4,
```

basis

```text
b(x) = (1, 2*x - 1, 6*x*(1-x)),
```

and strict-kernel amplitude **candidate**

```text
A_i(s,x) = sqrt(Cut_i(s)/ImM(s)) * K_strict(d_i(s,x)) * sqrt(x*(1-x)).
```

The exported tensor candidate is

```text
T_i,ab(s) = ∫_0^1 A_i(s,x)^2 b_a(x)b_b(x) dx.
```

The integral is evaluated by deterministic SciPy quadrature with exported error
matrices.

## Provenance boundary

The candidate is numerically executable, but it is not yet a strict-side
Feynman-derived backend amplitude.  Four nontrivial provenance gaps remain:

1. no exported strict Feynman-rule integrand,
2. no derivation of the loop measure from the strict action,
3. no derivation of the channel distance map from strict dynamics,
4. no normalization independent of the already exported `P2004` cut values.

Therefore the mathematical status of P2017 is

```text
OPEN_PROVENANCE_GAP_WITH_STRICT_QUADRATURE_TRACE
```

not a Cutkosky theorem closure.

## Machine-checkable result

The witness exports:

1. a tensor-candidate table for every `P2004` grid point and for every channel
   `gg/gh/hh/gx`,
2. quadrature error matrices,
3. eigenvalue checks for positive semidefiniteness,
4. trace positivity checks,
5. a coupled channel covariance candidate built from quadrature tensor traces,
6. a finite coupled scan using the quadrature-derived channel covariance
   candidate,
7. separate numerical and provenance gatekeeper blocks.

The numerical subdiagnostic may pass:

```text
PASS_STRICT_KERNEL_QUADRATURE_NUMERICS
```

but the exported theorem-level result remains:

```text
OPEN_PROVENANCE_GAP_WITH_STRICT_QUADRATURE_TRACE.
```

## No-false-pass boundary

P2017 is stronger than a pure tensor placeholder because it computes explicit
strict-kernel quadratures.  It is weaker than a backend theorem because its
amplitude is still a calibrated candidate rather than an independently derived
strict loop integrand.  It does not prove all-state unitarity, does not discharge
`QW-2191`, and does not close ToE.

## Next honest step

Build P2018 from explicit strict-side Feynman or loop-integrand rules: derive the
channel distance map, loop measure, normalization, and tensor integrand
independently of `P2004` cut calibration, then compare the resulting derived
tensors against the P2017 quadrature candidates on an expanded energy/RG atlas.
