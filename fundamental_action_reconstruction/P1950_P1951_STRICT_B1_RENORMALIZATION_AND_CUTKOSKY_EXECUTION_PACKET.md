# P1950/P1951 Strict B1 Renormalization And Cutkosky Execution Packet

Status: `OPEN_OBSTRUCTION_WITH_TRACE__LOCAL_B1_EXECUTION_ADVANCED_NO_FALSE_PASS`
As of: `2026-05-17`

## Goal

Execute the next strict-only numerical/symbolic step for:

1. B1 one-loop gravity divergence cancellation,
2. the first non-grid Cutkosky phase-space integral for
   `graviton->gauge_gauge`.

The packet uses only

```text
K_strict(d) = cos(omega*d + phi)/(1 + beta*d^eta)
```

with the QW-2049 tuple already used by `P1853`.

## P1950 Renormalization Result

`p1950_s900_strict_b1_renormalization_exact_cancellation.py` consumes the
declared `P1850/P1853` B1 coefficient layer and proves, in SymPy, that

```text
Gamma_div_B1 + delta_c_gr_1 R^2
             + delta_c_gr_2 Ricci^2
             + delta_c_gr_3 Riemann^2
             + delta_c_gr_4 G_GB = 0
```

after setting

```text
delta_c_gr_1 = -a_R2/epsilon
delta_c_gr_2 = -a_Ric2/epsilon
delta_c_gr_3 = -a_Riem2/epsilon
delta_c_gr_4 = -a_GB/epsilon
```

The local verdict is:

```text
PASS_ZERO_ON_DECLARED_B1_ANSATZ
```

This is a real algebraic cancellation theorem for the declared B1 ansatz.
It is not a first-principles derivation of the ansatz coefficients.

## P1951 Cutkosky Result

`p1951_s901_strict_b1_cutkosky_phase_space_integral.py` replaces the earlier
seed-grid proxy with the exact two-body massless phase-space integral for the
declared `P1860` seed integrand:

```text
K_cut_seed(s,x) = (a_R2 + a_Ric2*(1+x^2))*s,   x = cos(theta)
dPhi2 = dx/(32*pi)
```

Therefore

```text
Disc_seed(s) = s*(3*a_R2 + 4*a_Ric2)/(48*pi)
```

Since the evaluated strict B1 values have `a_R2 > 0`, `a_Ric2 > 0`, the seed
integral is nonnegative for `s > 0`.

The local verdict is:

```text
PASS_FOR_DECLARED_SEED_INTEGRAND_WITH_DOMAIN_OBSTRUCTION
```

## Hard Obstruction Found

`P1677` restricts the combined UR domain to:

```text
beta in [0.60, 0.90]
eta  in [1.10, 1.50]
```

But the current strict tuple used by `P1853` is:

```text
beta = 1.0
eta  = 1.8
```

So the `UR_link_theorem` cannot honestly be promoted from this packet. The
domain mismatch must be repaired by a new domain theorem, or the QW-2049
strict tuple must be handled by a separate UR object.

## Outputs

- `p1950_s900_strict_b1_renormalization_exact_cancellation.py`
- `generated/p1950_s900_strict_b1_renormalization_exact_cancellation_probe.json`
- `p1951_s901_strict_b1_cutkosky_phase_space_integral.py`
- `generated/p1951_s901_strict_b1_cutkosky_phase_space_integral_probe.json`

## No False Pass

This packet does not claim:

1. global renormalization closure,
2. full dressed Cutkosky equality,
3. full ghost-freedom,
4. background-independence globalization,
5. discharge of `QW-2191`.

It does claim that the B1 ansatz cancellation and the declared seed-integrand
phase-space integral are now executable, reproducible, and machine-checkable.

## Next Honest Step

Repair the `P1677` restricted-domain mismatch against QW-2049, then replace
the `K_cut_seed` integrand by a dressed amplitude whose imaginary part and cut
sum are reduced to the same symbolic basis.
