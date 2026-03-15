# F446 Current Strict Shannon Element‑Order Reference `O(2)`‑Cut Selector Ingredient Packet (No False‑PASS)

Status: `F446_EXECUTED_CURRENT_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_O2_CUT_SELECTOR_INGREDIENT_PACKET_NO_FALSE_PASS`  
As of: `2026-03-14`

## Goal

`T165` asks for an **actual** strict‑core Shannon‑typed selector ingredient that cuts the `QW-2191` continuous `O(2)`
family (on a declared `pair1` scope) down to an explicitly tracked residual discrete ambiguity, without smuggling a
marked direction/generator slot.

This packet exports one concrete ingredient instance:

```text
S_shannon_symmetry_breaking_selector_ingredient_o2_cut_v1
```

by pinning down:

1. the strict reference distribution `r_ord` built from `ord_Z12` + `alpha_geo_strict_derived_v1`,
2. the strict Shannon‑typed objective `J_ord` (cross‑entropy to `r_ord`),
3. the theorem‑level minimizer set (`θ = π/2 (mod π)`) from `N480`.

It does **not** claim strict per‑site vacuum/self‑coupling values (`T168`), does **not** decide `F2(d)` (`T166`), and
does **not** claim global `QW-2191` discharge or ToE closure.

## Strict‑admissible sources reused

1. `F329`
   - typed `Z_12` carrier language,
2. `F309/N420`
   - strict‑derived Shannon amplitude `alpha_geo_strict_derived_v1 := 4 ln 2`,
3. `N479`
   - `ord_Z12` is `Aut(Z_12)`‑invariant ⇒ no marked‑direction slot for `f(ord_Z12)` references,
4. `N480`
   - theorem‑level nontriviality + `Z2`‑unique minimizer for the `J_ord` objective on `pair1`.

## Exported ingredient (typed content)

### 1) Reference datum (element‑order weights; direction‑free)

Define:

$$
r_{\mathrm{ord}}(x)\propto \exp(-\alpha_{\mathrm{geo}}\,\operatorname{ord}_{Z_{12}}(x)),
\qquad
\alpha_{\mathrm{geo}}:=4\ln 2.
$$

Properties:

1. **No marked direction:** `r_ord` is `Aut(Z_12)`‑invariant (`N479`).
2. **Not translation‑invariant:** it marks the identity orbit `{0}` (explicit symmetry‑breaking premise class; see `T165`).

Persisted artifact:

`fundamental_action_reconstruction/generated/r_ord_z12_v1_reference_distribution.json`

### 2) Shannon‑typed selector objective (cross‑entropy)

On the `pair1` `O(2)` family (parameterized by $\theta$ as in `N480`), define:

$$
J_{\mathrm{ord}}(\theta):=-\sum_{x\in Z_{12}} p_\theta(x)\,\log r_{\mathrm{ord}}(x).
$$

By construction, `J_ord` depends on strict sources only through:

```text
alpha_geo_strict_derived_v1  +  (typed Z_12 -> ord_Z12)
```

and does not introduce any untracked generator/direction slot (`N479`).

### 3) Theorem‑level `O(2)` cut (residual `Z2`)

`N480` proves:

$$
\operatorname*{argmin}_{\theta\in[0,2\pi)} J_{\mathrm{ord}}(\theta)
=
\left\{\frac{\pi}{2},\frac{3\pi}{2}\right\}
=
\left\{\frac{\pi}{2}\ (\mathrm{mod}\ \pi)\right\}.
$$

So the induced cut is:

```text
O(2)  ->  Z2
```

with the residual `Z2` corresponding to the unavoidable sign flip `u_{θ+π} = -u_θ` (same squared‑amplitude distribution).

Persisted artifact:

`fundamental_action_reconstruction/generated/theta_fix_pair1_o2_cut_ord_reference_v1.json`

## Export name (as required by `T165`)

This packet instantiates the `T165` target object under the explicit reference datum:

```text
S_shannon_symmetry_breaking_selector_ingredient_o2_cut_v1
  := (r_ord, J_ord, theta_fix_set={pi/2 mod pi}, residual=Z2).
```

## What this does not claim

This packet does not claim:

1. strict per‑site `vpsi/g4/g6` provider (`T168`),
2. strict sigma six‑sum instantiation (`T167`),
3. strict diagonal defect decision (`T166`) or any strict diagonal/local `F2(d)` value,
4. strict-core selector closure or global `QW-2191` discharge,
5. ToE closure.

