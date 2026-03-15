# AX23 Strict‑Extension Lane Z12 Generator‑Fix + Shannon/KL O(2)‑Cut Premise Packet

Status: `AX23_EXECUTED_STRICT_EXTENSION_LANE_Z12_GENERATOR_FIX_AND_SHANNON_KL_O2_CUT_PREMISE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-14`

## Goal

`QW-2191` exports a strict-core continuous `O(2)` basis-choice family on degenerate 2D mode pairs (e.g. `pair1=span{c1,s1}`).

`T165` names the missing strict-core object-class for a Shannon-typed symmetry-breaking selector ingredient that would
canonically cut that family (up to an explicitly declared residual symmetry).

After `N463/N464`, any objective that depends only on squared site-amplitude distributions up to site permutation (in
particular any translation-invariant site-distribution functional) cannot yield a unique strict-core `O(2)` cut.

`P439` audits one explicit escape hatch: a **non-translation-invariant** KL-to-reference objective using a marked-site +
marked-direction reference distribution on the 12-site ring. It finds a **Z2-unique** minimizer pattern (two minima
separated by `π`) on a dense grid for one directed reference `r_dir(x) ∝ exp(-(4 ln 2)*x)`.

However, `N462` proves there is no `Aut(Z_12)`-invariant canonical generator/orientation selection in strict core from
typed `Z_12` + `Aut(Z_12)` structure alone. Therefore any directed reference `r_dir(x)` implicitly fixes a generator
(`x -> x+1`) and cannot be promoted as strict-core canonical without discharging `T164` (or restricting the admissible
symmetry in a tracked way).

So the next honest move **today** (if one insists on a single reproducible representative of the `O(2)` family) is not
strict-core promotion. It is an explicit symmetry-breaking **premise** in the already accepted separated scope:

```text
strict_extension_only
```

`AX23` records that premise explicitly and freezes:

1. a `Z_12` generator/orientation choice, and
2. one Shannon/KL-based `pair1` representative angle `theta_pair1`,

without changing strict core and without implying `T165` / `QW-2191` discharge.

## Inputs reused

1. `AX16`
   - `strict_extension_only` scope is accepted (theory-level).
2. `F329/N450` + `F331/N453`
   - typed `Z_12_v1` carrier and `Aut_Z12_v1` ambiguity group.
3. `N462` + `T164`
   - strict-core nonexistence of an Aut-invariant canonical generator/orientation choice, and the named missing fixing-datum target.
4. `N420`
   - strict Shannon source upgrade:
     `alpha_geo_strict_derived_v1 := 4 ln 2`.
5. `P439`
   - probe-level audit that the directed reference `r_dir(x) ∝ exp(-(4 ln 2)*x)` yields a Z2-unique minimizer pattern
     for `KL(u1(theta)^2 || r_dir)` on a dense grid.

## What is accepted (extension scope only)

`AX23` accepts the following statement **only** in `strict_extension_only` scope:

1. **Fix a preferred generator/orientation of `Z_12_v1`:**
   - choose the successor step `suc_fix(k) := (k+1) mod 12` (equivalently: generator `g_fix := 1 ∈ Z_12_v1`),
   - and read the directed coordinate of a site as the integer label `x ∈ {0,..,11}` relative to the identity.

2. **Define a directed Shannon/KL reference distribution on ring sites:**
   $$
   r_{\mathrm{dir}}(x)
   := \frac{\exp(-\alpha_{\mathrm{geo}}\,x)}{\sum_{t=0}^{11}\exp(-\alpha_{\mathrm{geo}}\,t)},
   \qquad
   \alpha_{\mathrm{geo}} := 4\ln 2.
   $$

3. **Define the `pair1` objective on the `QW-2191` family:**
   - `u1(theta) := cos(theta)c1 + sin(theta)s1`,
   - `p_theta(x) := u1(theta)(x)^2`,
   - objective: `J_dir(theta) := KL(p_theta || r_dir)`.

4. **Freeze one reproducible `pair1` representative angle:**
   Choose `theta_pair1_fix := π/6` as the representative in `[0,π)`, justified as the Z2-unique minimizer pattern
   observed by `P439` (two minima separated by `π`).

This is a *premise-based selection rule*, not a strict-core derivation, not a uniqueness theorem, and not a discharge of `T165`.

## Meaning

This packet means only:

1. the repo can now cite one reproducible `pair1` representative angle from a Shannon/KL selector **in extension scope**,
2. strict core remains unchanged and honest: `T165` stays open and `QW-2191` stays undischarged in strict core,
3. any future strict-core promotion must discharge the missing fixing-datum target (`T164`) and the typed selector ingredient target (`T165`).

## Hard limits (no false pass)

`AX23` does not claim:

1. strict-core canonical generator/orientation fixing datum export (`T164` remains open),
2. strict-core Shannon selector ingredient export (`T165` remains open),
3. strict-core theta export,
4. strict-core selector closure or strict-core `QW-2191` discharge,
5. ToE closure.

