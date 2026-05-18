# P1969 / S919 — Strict B1 Exact Phase-Space Integral First Non-Proxy Theorem Packet

Status: `OPEN_OBSTRUCTION_WITH_TRACE`  
Route: `strict_only`  
Date anchor: `2026-05-18`

## Goal

Execute the first honest non-grid Cutkosky calculation after `P1864` by replacing
the pointwise `P1860` projected discontinuity samples with an exact two-body
phase-space integral for the same declared projected B1 polynomial backend.

This packet is intentionally narrow. It does **not** claim full unitarity,
full dressed-residue closure, RG-stable pole transport, background-family
closure, or ToE closure.

## Strict inputs

The strict coefficient source is `P1853`:

```text
K_strict(d) = cos(omega*d + phi)/(1 + beta*d^eta)
omega = 0.18575, phi = 0.16250, beta = 1.0, eta = 1.8
alpha_geo_strict = 4 ln(2)
```

The two coefficients used from the B1 layer are:

```text
a_R2   = 18516431*log(2)/(640000000*pi**2)
a_Ric2 = 9937*log(2)/(40000*pi**2)
```

No legacy-kernel role transfer is used.

## Integrated object

`P1860` declared the projected polynomial sample

```text
K_cut_sample(s, theta) := (a_R2 + a_Ric2*(1 + theta^2))*s.
```

`P1969` promotes the angular variable to

```text
x = cos(theta),  x in [-1, 1],
```

and integrates

```text
K_proj(s,x) = s*(a_R2 + a_Ric2*(1+x^2))
```

with the massless two-body convention

```text
dPi_2 = (1/(32*pi^2)) dphi dx,
phi in [0, 2*pi], x in [-1, 1].
```

## Symbolic result

Sympy evaluates

```text
Disc_proj_B1(s)
= ∫ dPi_2 K_proj(s,x)
= s*(a_R2/(8*pi) + a_Ric2/(6*pi)).
```

The exported residual check is exactly zero:

```text
Disc_proj_B1(s) - s*(a_R2/(8*pi) + a_Ric2/(6*pi)) = 0.
```

Because `a_R2 > 0`, `a_Ric2 > 0`, and `s > 0`, this projected B1 integral is
positive in the declared physical positive-`s` corridor.

## Machine-check export

The execution script is:

```text
fundamental_action_reconstruction/p1969_s919_strict_b1_exact_phase_space_integral_first_nonproxy_checkpoint.py
```

It writes:

```text
fundamental_action_reconstruction/generated/p1969_s919_strict_b1_exact_phase_space_integral_first_nonproxy_checkpoint.json
```

The JSON contains:

1. the symbolic closed form,
2. the zero residual check,
3. an independent `scipy.integrate.dblquad` cross-check on `s = 0.5, 1, 2, 4, 8`,
4. conservative integration-error lower bounds,
5. explicit no-false-pass boundaries.

## Theorem-level statement actually licensed

Within the `P1860` projected B1 polynomial backend and the stated massless
two-body phase-space convention,

```text
Disc_proj_B1(s) = s*(a_R2/(8*pi) + a_Ric2/(6*pi)) > 0
```

for `s > 0` under the `P1853` strict coefficients.

This is a real symbolic/numeric improvement over a finite point-sample table:
it integrates the declared projected polynomial over the full two-body angular
phase space.

## Honest boundary

This is **not** the requested full `UR_link_theorem` yet. The following remain
open:

1. full dressed graviton propagator pole computation,
2. full `graviton -> gauge_gauge` numerator and polarization sum,
3. positive-residue certification for the full dressed amplitude,
4. RG transport of allowed poles and residues,
5. background-family continuation beyond B1,
6. global unitarity, TG3 closure, and ToE closure.

Therefore the status remains:

```text
OPEN_OBSTRUCTION_WITH_TRACE
```

## Next honest step

Compute the dressed `graviton -> gauge_gauge` amplitude numerator and pole
residues, then rerun this exact phase-space pipeline on the full
polarization-summed integrand instead of the `P1860` projected polynomial.

## Lay explanation

Dotąd teoria miała kilka przykładowych punktów testu unitarności. Teraz dla
jednego uproszczonego, ale jawnie zadeklarowanego jądra policzyliśmy całą całkę
po kierunkach dwóch cząstek końcowych. To jest postęp, bo przechodzimy od
„tabelki punktów” do prawdziwej całki. Nie jest to jednak jeszcze pełny dowód:
brakuje pełnej amplitudy fizycznej i dokładnych biegunów propagatora.
