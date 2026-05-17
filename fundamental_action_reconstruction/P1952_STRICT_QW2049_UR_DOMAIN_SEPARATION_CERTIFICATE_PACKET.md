# P1952 Strict QW-2049 UR Domain Separation Certificate Packet

Status: `OPEN_OBSTRUCTION_WITH_TRACE__QW2049_SEED_LOCAL_DOMAIN_CERTIFIED_NO_FALSE_PASS`
As of: `2026-05-17`

## Goal

Respond to the obstruction found in `P1951`:

```text
P1677 beta/eta domain does not contain the QW-2049 strict tuple.
```

The correct repair is not to silently widen `P1677`, but to export a separate
QW-2049 seed-local domain certificate.

## Domain Result

`P1952` exports the QW-2049 seed-local rectangle:

```text
beta  in [19/20, 21/20]
eta   in [17/10, 19/10]
omega in [9/50, 19/100]
phi   in [3/20, 17/100]
```

This rectangle contains the QW-2049 strict tuple:

```text
beta = 1
eta = 9/5
omega = 743/4000
phi = 13/80
```

It is disjoint from the old `P1677` domain:

```text
beta in [0.60, 0.90]
eta  in [1.10, 1.50]
```

## Positivity Result

On the QW-2049 seed-local rectangle, `P1952` proves strict positive lower
bounds for:

```text
a_R2
a_Ric2
a_Riem2
a_GB
Disc_seed(s)/s
```

Therefore, for `s > 0`,

```text
Disc_seed(s) > 0
```

on that seed-local rectangle.

## Safe Statement

The declared B1 counterterm cancellation from `P1950` and the seed Cutkosky
positivity from `P1951` coexist on a certified QW-2049 seed-local domain.

## Rejected False Pass

`P1952` does not claim:

1. extension of `P1677`,
2. full `UR_link_theorem`,
3. dressed graviton-to-gauge-gauge Cutkosky equality,
4. background-global renormalization,
5. BRST physical-state projection,
6. `QW-2191` discharge.

## Outputs

- `p1952_s902_strict_qw2049_ur_domain_separation_certificate.py`
- `generated/p1952_s902_strict_qw2049_ur_domain_separation_certificate.json`

## Next Honest Step

Build `P1953`: replace the seed integrand by a dressed amplitude model reduced
to the same symbolic basis, or export a precise obstruction if the dressed
amplitude is not available in the current repo.
