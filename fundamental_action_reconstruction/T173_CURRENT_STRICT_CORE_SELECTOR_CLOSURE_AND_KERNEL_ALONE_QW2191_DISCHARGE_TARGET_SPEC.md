# T173 Current Strict-Core Selector Closure + Kernel-Alone Global QW-2191 Discharge Target Spec

Status: `T173_CURRENT_STRICT_CORE_SELECTOR_CLOSURE_AND_KERNEL_ALONE_QW2191_DISCHARGE_TARGET_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

After the repo now exports:

1. global selector atlas + transition/gluing objects on `C_v1` (`T170` discharged: `F469/N515`),
2. global selector state objects on `C_v1`:
   - projective (`F470/N516`),
   - directed in an explicit premise-based scope (`T171` discharged: `F473/N523`, `F474/N524`),
3. a fully admissible strict-core source object for the internal selector source step `S_sel_int` (full `F34` contract discharged: `N676`),
4. global seed‑v1 operator-chain promotions on `C_v1` (`N550–N553`),
5. global selector closure objects on `C_v1` in both scopes (`T172` discharged in:
   - projective scope: `N674` (via `F672/N672/N673`),
   - directed scope: `N678` (via `F677/N677`), with the raw-output obstruction boundary recorded by `N675`),
6. a theorem-level **projective strict-core selector closure** discharge (`N680`),
7. and a theorem-level boundary that packages what remains beyond the closure objects (`N679`),

the remaining honest strict frontier is now **post‑`T172`**:

```text
kernel-alone/global QW-2191 discharge claim discipline
and any directed/sign-sensitive physical orientation datum in strict core beyond projective/ray semantics
```

`T173` exists to prevent a false conflation:

- `T172` discharges **closure objects** (projective/directed) and their well-definedness in declared scopes,
- `N680` discharges **projective strict-core selector closure** (ray-level) post-`T172`,
- but `T173` remains the target label under which one may even *attempt* to claim:
  - `QW2191_kernel_alone_discharge = true`, and/or
  - any directed/sign-sensitive physical orientation datum in strict core without smuggled premises.

This is a **target spec only**. It exports no new object.

## Scope

`T173` is scoped to the strict-core selector-closure frontier and kernel-alone/global `QW-2191` discharge discipline.

`T173` does **not** decide:

1. ToE closure,
2. any legacy→strict role transfer beyond what is explicitly exported,
3. any operator-level transition groupoid promotion beyond projector/section level (`N512` boundary).

## Target objects (what would count as a discharge)

### A) Strict-core selector closure (closure-level theorem, declared scope)

Export one theorem-level package (and its audited generated summary) that pins down, explicitly:

1. what is meant by “strict-core selector closure” in the declared scope (projective vs directed),
2. that the closure outcome is uniquely determined **from exported strict-core objects** in that scope (no hidden selector knobs),
3. which premises are tracked (e.g. `T164` generator/orientation fixing datum; explicit sign-lift conventions if used),
4. that the result does **not** smuggle `QW-2192` (selector axiom) into strict core.

Minimal acceptable form:

```text
StrictCoreSelectorClosure_global_C_v1_strict_v1 (theorem-level package):
  - strict_core_selector_closure = true
  - QW2191_kernel_alone_discharge = false (unless separately proven)
  - ToE_closure = false
```

Current-state note: projective strict-core selector closure is now discharged by `N680` with scope label
`projective_ray_state` (see `fundamental_action_reconstruction/generated/n680_current_strict_t173_projective_strict_core_selector_closure_discharge_theorem_summary.json`).

### B) Kernel-alone/global `QW-2191` discharge or boundary (theorem-level)

Export one explicit theorem-level statement that makes precise what is meant by:

- “kernel-alone/global `QW-2191` discharge”, and
- whether it is:
  1) actually discharged (uniqueness from kernel alone in the declared strict scope), or
  2) impossible / still open as a strict boundary.

Any positive claim here must keep the distinction explicit:

```text
QW-2191 (kernel-alone obstruction) is not "false" unless a genuine theorem-level uniqueness statement is exported.
```

## Acceptance tests (no false pass)

Any honest `T173`-class discharge must satisfy all of:

1. **Explicit scope:** projective vs directed must be declared; any directed claim must keep sign-lift / fixing-datum dependence explicit (`N462/T164`, `N675`).
2. **Global well-definedness:** must rely only on exported global atlas/transition/state/operator objects on `C_v1` (no “implicit overlap” semantics).
3. **Level discipline:** any cocycle/path-independence claim must respect the projector/section-level boundary (`N512`).
4. **No selector smuggling:** no silent import of `QW-2192` into strict core.
5. **No implied ToE closure:** ToE closure remains separate and must stay false unless separately discharged.

## Current-state note

On the current repo state:

- `T172` is discharged at the closure-object level (`N674`, `N678`),
- projective strict-core selector closure is discharged at theorem level (`N680`),
- while `QW2191_kernel_alone_discharge = false` and any directed/sign-sensitive physical orientation datum in strict core remain explicitly unclaimed,
  packaged by the frontier boundary theorem `N679`.

Update (`2026-03-17`): even when attempting to construct a **deterministic** per-chart sign-lift rule from the exported strict-core
seed payload weights (`F647`), the directed output sign is not globally chart-independent (`P681`), packaged as a strict boundary by `N681`.
A chart-sine-aligned sign-lift **convention** can be made output-sign-consistent (`P682`), but this depends on a non-`Aut(Z_12)`-invariant chart embedding
and therefore does **not** upgrade any directed/sign-sensitive physical orientation datum into strict core.
Moreover, a rooted transport-based sign lift can be made output-sign-consistent by fixing sign on `pair1` from the exported reflection-breaking weight `w_break`
and propagating to other charts via the exported `O_1m` transports (`P683`), but this still depends on axis-only transport representatives and remains a section/convention choice.
This rooted sign-lift is also exported as an explicit global sign-lift/section-choice object on `C_v1` (`F683`), without any physical orientation claim.
Finally, the same `w_break`-rooted rule is now also exported as a global **directed state representative** in an explicit convention scope on `C_v1`
(`F684`, audited by `P684`), descending to the strict projective state but still below any physical sign datum claim.

Therefore, `T173` is the next honest strict label for continued closure attempts beyond `T172`.

## Hard limits

`T173` must not claim:

1. ToE closure,
2. global kernel-alone `QW-2191` discharge unless separately proven,
3. any `Aut(Z_12)`-invariant sign canonicity “for free”.
