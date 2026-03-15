# P446 Current Strict `T168` Value‑Provider Real‑Discharge Blueprint After `F446` (`T165`) (No False‑PASS)

Status: `P446_EXECUTED_CURRENT_STRICT_T168_VALUE_PROVIDER_REAL_DISCHARGE_BLUEPRINT_AFTER_F446_T165_NO_FALSE_PASS`  
As of: `2026-03-14`

## Goal

After `P444/P438`, the repo is explicit:

```text
T168 is the real strict bottleneck under the diagonal/local accelerator lane.
```

After `F446` + `N480`, the repo now also has one **direction‑free** strict Shannon‑typed selector ingredient:

```text
r_ord(x) ∝ exp(-alpha_geo*ord_Z12(x))  +  J_ord = cross_entropy(p_theta || r_ord)
  ->  θ = π/2 (mod π)  (pair1 O(2) -> Z2 cut).
```

This pointer records the most honest *next strict move* for actually discharging `T168` **without** smuggling a
marked‑direction slot — i.e. using `F446` only as an explicit ingredient, not as a verbal story.

It does **not** discharge `T168`.

## Why `T168` still does not follow from current strict scalar closure

`N478` already proves: strict scalar vacuum closure (`QW-2122/2123/2124`) does not canonically lift to per‑site arrays
`vpsi[0..11], g4[0..11], g6[0..11]`.

So `T168` needs an additional strict ingredient that supplies a **canonical lift** (mapping/selector) from strict scalar
data to **per‑site** data.

`F446` supplies a clean, direction‑free *candidate class* of such an ingredient, but by itself it does not yet define a
per‑site vacuum/self‑coupling provider.

## The only realistic strict pattern: “constrained uniqueness”

To make a `T168` provider genuinely `strict_derived`, the repo must export (as strict objects) a chain of the form:

```text
(strict scalar vacuum data) + (strict kernel channel) + (strict selector ingredient)
  -> unique per-site arrays (vpsi,g4,g6)
```

Concretely, the next strict work must produce **all** of the following, with no hidden slots:

1. a **declared candidate family** of admissible per‑site arrays (domain),
2. a **strict objective / extremality principle** that selects a unique member of that family,
3. a **theorem‑level existence + uniqueness** proof (or a strictly bounded discrete residual ambiguity),
4. an exported **value object** (JSON) containing the resulting numeric arrays.

Anything less remains either:

- underdetermined (and cannot be called `strict_derived`), or
- premise‑based (must remain `strict_extension_only`).

## A concrete strict‑compatible blueprint (uses `F446` explicitly, no direction slot)

This is the cleanest blueprint currently compatible with the repo’s “no false‑PASS” discipline:

### Step 1 — Fix strict inputs already exported

Reuse strict objects already present:

1. `rho_star_sq`, `lambda_psi_strict`, branch rule (from `QW-2122/2123/2124`),
2. strict kernel channel `K_total` and diagonal floor `m0^2` (`R14/R15`),
3. the direction‑free reference datum `r_ord` from `F446`.

### Step 2 — Declare an admissible per‑site candidate domain (no hidden origin/direction)

Define a candidate set `C_vac` of per‑site arrays.

At minimum it must:

1. encode the scalar norm constraint (e.g. `Σ_i vpsi_i^2 = rho_star_sq`),
2. make explicit whether the nonzero premise `vpsi_i ≠ 0` is required (for `N474/N475/N477`),
3. state all symmetry constraints you intend to keep (e.g. `Aut(Z_12)`‑invariance is admissible; translation‑invariance
   is **not** if you want `F2(d)≠0`),
4. forbid “choose a site / choose a generator / choose a representative” slots.

### Step 3 — Export a strict selector functional on that domain

Define an explicit functional `J_vac` on `C_vac`.

The minimal direction‑free Shannon‑typed candidate is the **element‑order expectation** induced by `F446`:

$$
J_{\mathrm{ord}}(v):=\sum_{x\in Z_{12}} q_v(x)\,\operatorname{ord}(x),
\qquad
q_v(x):=\frac{v(x)^2}{\sum_y v(y)^2}.
$$

But **alone** this will select a trivial extreme (concentrate on the identity orbit), so it must be combined with the
actual canonical vacuum constraints (stationarity) rather than used as a free unconstrained fit.

Therefore the strict version must be stated as:

```text
minimize J_ord(v) subject to: (canonical constant-vacuum stationarity constraints) + (scalar closure constraints).
```

If this is the chosen route, the repo must then export:

1. the exact constraint system,
2. an existence proof,
3. a uniqueness proof (or explicit residual finite ambiguity),
4. the resulting numeric arrays as a value object consumable by `P437`.

### Step 4 — Produce the strict value provider object

Export one JSON value object classified as `strict_derived` that provides:

```text
vpsi[0..11], g4[0..11], g6[0..11]
```

and explicitly references its derivation chain (constraints + selector functional + uniqueness theorem).

Only then `P437 -> T167 -> P434 -> T166` becomes strictly computable.

## What to do next in practice (repo work items)

The next honest work items under this blueprint are:

1. author a new strict target spec that makes the above “constrained uniqueness” mapping explicit (either extend `T168`
   or add a child target). (`T169` now names this lift explicitly.) and
2. implement a probe‑level solver harness (explicitly labeled) to test whether the constrained problem even has a
   nontrivial solution on current exports, **before** attempting theorem‑level promotion.

## What P446 does not claim

`P446` does not claim:

1. that the above constrained uniqueness mapping is already exported,
2. that `T168` is already discharged,
3. that `F2(d)≠0` holds in strict core,
4. global `QW-2191` discharge,
5. ToE closure.
