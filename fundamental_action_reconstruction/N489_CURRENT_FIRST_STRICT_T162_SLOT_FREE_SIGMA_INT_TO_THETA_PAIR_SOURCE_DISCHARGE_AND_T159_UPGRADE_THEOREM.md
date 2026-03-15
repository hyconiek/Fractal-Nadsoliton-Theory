# N489 Current First Strict `T162` Slot‑Free Sigma‑Int → Theta‑Pair Source Discharge and `T159` Upgrade Theorem

Status: `N489_DISCHARGED_CURRENT_FIRST_STRICT_T162_SLOT_FREE_SIGMA_INT_TO_THETA_PAIR_SOURCE_DISCHARGE_AND_T159_UPGRADE_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

`T159` names a strict‑core upgrade target: export a canonical sigma‑int → theta selector ingredient whose `O(2)` cut is
strict‑core (no hidden selector slots).

On the current repo state, the old sigma‑int → theta candidate class has two exposed selector slots (`eps`, `delta_d`)
and cannot be upgraded by invariance (`N443`). Strict‑derived slot selection (`T160/T161`) is not exported.

Therefore the only honest strict‑core continuation is the construction‑class change route `T162`.

This theorem records the first strict discharge of that route:

```text
export a slot‑free sigma_int -> (theta_1, theta_2) theta‑pair source,
derived from the strict Shannon element‑order reference objective class,
thereby discharging T162 and satisfying T159 acceptance test (2C).
```

## Strict‑admissible evidence reused

1. `T159`
   - strict‑core upgrade target spec and acceptance tests,
2. `T162` + `N449`
   - slot‑free construction class target is named and was not previously discharged,
3. `F446` + `N480`
   - strict, direction‑free Shannon selector ingredient on `pair1`:
     `argmin J_ord,1(θ) = {π/2 (mod π)}` (`O(2)->Z2` cut),
4. `N488`
   - strict theorem‑level extension to `pair2`:
     `argmin J_ord,2(θ) = {0 (mod π)}` (`O(2)->Z2` cut),
5. `F451`
   - actual exported slot‑free sigma‑int → theta‑pair source object:
     `ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1`.

## Theorem (discharge of `T162` and satisfaction of `T159` via 2C)

### Claim 1. An actual slot‑free sigma‑int → theta‑pair source object is exported.

`F451` exports an actual strict theta‑pair source artifact:

```text
ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1
  : sigma_int_strict_derived_v1 -> (theta_1, theta_2)
```

with:

1. no `eps` amplitude parameter,
2. no `delta_d` corridor‑step parameter,
3. no theta inputs and no populated basis‑pair inputs (noncyclic contract),
4. observer‑free (no `K_obs` indexing).

So a **slot‑free construction class instance exists** in strict core. ∎

### Claim 2. The exported theta‑pair cuts `pair1` and `pair2` from continuous `O(2)` to residual `Z2`.

By `N480`, on `pair1` the strict Shannon element‑order reference cross‑entropy objective has minimizer set:

$$
\arg\min_{\theta\in[0,2\pi)} J_{\mathrm{ord},1}(\theta)=\left\{\frac{\pi}{2}\ (\mathrm{mod}\ \pi)\right\}.
$$

By `N488`, on `pair2` the same objective class has minimizer set:

$$
\arg\min_{\theta\in[0,2\pi)} J_{\mathrm{ord},2}(\theta)=\left\{0\ (\mathrm{mod}\ \pi)\right\}.
$$

Each minimizer set is exactly residual `Z2` (sign flip `u_{θ+π}=-u_θ` yields the same squared‑amplitude distribution).

`F451` fixes explicit representatives `(theta_1,theta_2)` consistent with those theorem‑level sets, and it explicitly
tracks the residual `Z2` on `pair1` via the already exported sigma‑int sign convention (`F311`).

Therefore, in the declared scope (the `R1` residual orientation datum slot which depends on `pair1` and `pair2`), the
continuous `O(2)` ambiguity is reduced to an explicitly tracked residual discrete ambiguity. ∎

### Claim 3. `T162` is discharged.

`T162` asks for an exported sigma‑int → theta construction class with **no** `eps`/`delta_d` selector slot families.

By Claim 1, such a slot‑free object is exported (`F451`), with explicit theorem‑level uniqueness/cut arguments on the
relevant pair planes (`N480`, `N488`).

Therefore the `T162` slot‑free construction‑class target is discharged on the current repo state. ∎

### Claim 4. `T159` is satisfied via acceptance test (2C).

`T159` acceptance test (2C) allows a strict‑core upgrade by changing the construction class so the selector slots do not
exist at all, provided the result supplies a canonical `O(2)` cut in strict core (up to residual `Z2`) with no hidden
slots.

By Claim 3, the slot‑free class exists and is exported.
By Claim 2, it supplies a theorem‑level `O(2)->Z2` cut on `pair1` and `pair2` (the `R1` scope).

Therefore the strict sigma‑int → theta strict‑core upgrade target `T159` is satisfied in the stated scope by the
discharged slot‑free construction class. ∎

## What this theorem does not claim

`N489` does not claim:

1. strict‑core selector closure (`S_sel_int`) or admissibility,
2. global `QW-2191` discharge beyond the stated `pair1/pair2` scope,
3. discharge of the eps/delta slot‑selection targets (`T160/T161`),
4. discharge of post‑`T148` object‑support targets (`N302/N395/T130`),
5. ToE closure.

