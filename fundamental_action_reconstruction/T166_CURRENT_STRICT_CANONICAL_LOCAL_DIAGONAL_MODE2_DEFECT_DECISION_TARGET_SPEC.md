# T166 Current Strict Canonical Local Diagonal Mode‑2 Defect Decision Target Spec

Status: `T166_CURRENT_STRICT_CANONICAL_LOCAL_DIAGONAL_MODE2_DEFECT_DECISION_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-13`

## Goal

After `QW-2191`, one concrete strict candidate “accelerator of choice” on `pair1 = span{c1,s1}` is a **diagonal/local**
sector breaking the continuous `O(2)` basis-choice family.

`N466` reduces any diagonal/local `pair1` `O(2)` cut attempt to one explicit mode‑2 defect:

$$
D=\mathrm{diag}(d_0,\ldots,d_{11})\ \text{cuts `O(2)` on `pair1`}\quad\Longleftrightarrow\quad
F_2(d)\neq 0,
$$

where (for `n=12`):

$$
F_2(d):=\sum_{k=0}^{11} d_k\,e^{i\frac{4\pi k}{12}}\in\mathbb{C}.
$$

`N467/P426` then reduce the **canonical FIN local diagonal residual sector** question to an explicit six‑class linear
combination.

`T166` names the missing strict object:

```text
strict-derived decision of F2(d) for the canonical FIN D_local_residual
  (either prove F2(d)=0, or prove F2(d)≠0),
without hidden site/generator slots.
```

This target is intentionally narrow: it is the minimal strict step that would either (a) kill the diagonal accelerator
route on `pair1`, or (b) turn it into a checkable strict `pair1` `O(2)`-cut ingredient.

## Scope

`T166` is scoped only to the canonical FIN diagonal/local residual sector on the strict `n=12` carrier and to the single
mode `pair1`.

It does **not** decide:

1. any strict-core theta export (`T159`) beyond the `pair1` axis,
2. any strict eps/delta slot elimination (`T160/T161`),
3. any slot-free sigma-int → theta class (`T162`),
4. strict-core selector closure or global `QW-2191` discharge,
5. ToE closure.

## Target output objects

Discharging `T166` must export at least one of:

1. a strict theorem:
   ```text
   N*_canonical_local_diagonal_mode2_defect_zero :  F2(d)=0
   ```
   (diagonal/local `pair1` accelerator route is dead), **or**
2. a strict theorem:
   ```text
   N*_canonical_local_diagonal_mode2_defect_nonzero :  F2(d)≠0
   ```
   (diagonal/local `pair1` accelerator route is live), together with the induced canonical `pair1` axis datum
   computable by `N468`.

## Strict inputs required (must be explicit)

Any discharge must explicitly reuse:

1. `R15`
   - the canonical diagonal decomposition:
     ```text
     D_canonical = m0^2 I + D_local_residual
     ```
2. `N466`
   - diagonal `pair1` `O(2)`-cut criterion via `F2(d)`,
3. `N467/P426`
   - the six‑class reduction of `F2(d)` for `n=12` in terms of the opposite‑pair sums
     $S_k:=d_k+d_{k+6}$ (or the exported `R18` class language).

If the discharge uses additional structure to fix the diagonal profile (e.g. vacuum/EOM solving, micro constraints,
route-scoped identifications), it must name those inputs explicitly and prove they are strict‑admissible in the declared
scope.

## Acceptance tests (what counts as discharge)

An actual discharge of `T166` must provide:

1. **Canonical profile identification:** a strict statement pinning down which diagonal profile $d_k$ is being decided,
   namely the canonical FIN residual profile:
   $d_k = (D_{\mathrm{local\_residual}})_{kk}$ from `R15`.
2. **No hidden slot:** the discharge must not introduce an untracked “choose a site”, “choose a generator”, or “choose a
   gauge representative” slot to decide $F_2(d)$.
3. **Strict-derived constraints or values:** export strict-derived relations sufficient to decide the two reduced
   constraints:
   $$
   \mathrm{Re}\,F_2 = 0,\qquad \mathrm{Im}\,F_2 = 0
   $$
   (for the `=0` case), or prove at least one is nonzero (for the `≠0` case), *for the canonical profile*.
4. **Nontriviality guard:** it is **not** sufficient to exhibit an *ad hoc* assignment inside the free `R15` coefficient
   class; such moves are closed as underdetermined by `N472/P431`.
5. **If `F2(d)≠0`: canonical axis exportability:** the discharge must explicitly invoke `N468` (or re-prove its
   conclusion) to state how the unique `pair1` axis/eigenbasis follows (up to the declared residual `Z2`).
6. **No false pass discipline:** the discharge must not claim strict-core selector closure, global `QW-2191` discharge,
   or ToE closure unless separately proved.

## Current status (why `T166` is still open)

On the current repo state, `N472/P431` prove that at the level of the exported `R15` canonical diagonal coefficient
class, $F_2(d)$ is **underdetermined**: both $F_2(d)=0$ and $F_2(d)\neq 0$ remain compatible with current exports.

Therefore `T166` is not discharged today and must be attacked by exporting additional strict-derived structure that
decides the canonical diagonal residual profile.

Additionally, `P432` audits that no already-exported **numeric value instantiation** exists in the repo that would
decide the `P426` expression anyway. So `T166` remains open both:

1. at the coefficient-class level (`N472/P431` underdetermination), and
2. at the repo-state “value instantiation present?” level (`P432` scan).
