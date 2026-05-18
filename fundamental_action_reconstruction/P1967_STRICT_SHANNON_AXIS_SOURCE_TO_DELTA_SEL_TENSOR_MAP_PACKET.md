# P1967 Strict Shannon Axis Source to `Delta_sel` Tensor Map Packet

Status: `STRICT_SHANNON_AXIS_SOURCE_INSTANTIATES_P1966_DELTA_SEL_ON_PAIR1_TO_PAIR5__GLOBAL_SELECTOR_CLOSURE_STILL_OPEN`  
As of: `2026-05-17`

## Goal

Address the repository-grep correction after `P1966`.

`P1966` correctly identified the minimal tensor shape required to cut a
Fourier-pair `O(2)` family, but its initial audit was too coarse: the repo
already contains strict, lane-scoped **axis-only** selector sources.  The next
honest move is therefore not to keep saying “no source exists”, but to map the
existing strict source into the exact `P1966` tensor form and preserve the
remaining limits.

## Grep result used

The relevant strict-lane source family found in the repo is:

1. `N480` — Shannon element-order reference cuts `pair1` from `O(2)` to
   residual `Z2`.
2. `N488` — the same source cuts `pair2`.
3. `N496` — the same source cuts `pair3..pair5`.
4. `F454` — exports the mode-index assignment object
   `ModeIndexAssignment_shannon_element_order_reference_strict_core_v1`.
5. `N514` — packages nonvanishing of the element-order Fourier defect.
6. `P455` — audits alignment between the Shannon source and the diagonal/local
   source up to residual `Z2`.
7. `B8` — keeps the no-false-pass boundary: no global `QW-2191` discharge and
   no admissible `S_sel_int` closure.

## Construction

On the strict `n=12` scaffold define

```text
ord_Z12(0) = 1,
ord_Z12(x) = 12/gcd(x,12) for x != 0.
```

For each degenerate pair `pair_m`, `m=1..5`, compute

```text
F_2m(ord) = sum_{x=0}^{11} ord_Z12(x) * exp(i*4*pi*m*x/12).
```

The `P1966` minimal tensor is then instantiated by

```text
Delta_sel_pair_m = [[Re(F_2m), Im(F_2m)],
                    [Im(F_2m), -Re(F_2m)]].
```

The spectral gap is

```text
gap_m = 2*abs(F_2m).
```

Therefore `F_2m != 0` is exactly the condition that the `O(2)` family is cut
to residual `Z2` on `pair_m`.

## Machine result

The executable probe `p1967_s917_strict_shannon_axis_source_to_delta_sel_tensor_map.py`
computes the five defects directly and compares them to the exported `F454`
mode-index assignment witness.  The generated JSON records:

```text
all_pair_gaps_nonzero = true
all_pair_defects_match_exported_F454 = true
p455_shannon_vs_diagonal_axis_alignment_up_to_z2 = true
strict_axis_only_source_pass = true
```

Thus the existing strict Shannon element-order lane does instantiate the
minimal `P1966` symmetry-breaking tensor on `pair1..pair5`.

## No false pass

This packet proves only a lane-scoped axis result:

```text
O(2) -> residual Z2 on the exported n=12 pair planes.
```

It does **not** prove:

1. sign-sensitive physical orientation,
2. admissible `S_sel_int`,
3. global `QW-2191` discharge,
4. `F(nadsoliton) => L_SM + L_GR` mapping witness,
5. full ToE closure.

## Output artifacts

- Script:
  `p1967_s917_strict_shannon_axis_source_to_delta_sel_tensor_map.py`
- Witness:
  `generated/p1967_s917_strict_shannon_axis_source_to_delta_sel_tensor_map.json`

## Next honest step

Test whether the downstream observables needed for strict ToE closure are
insensitive to the residual `Z2` sign on the `P1967` axes.  If they are, export
a gauge-irrelevance theorem for those observables.  If not, construct or refute
a strict sign-sensitive datum.

## Lay explanation

The repo already has a strict way to choose an **axis** in each ambiguous
2D plane: use the Shannon element-order profile on `Z12`.  `P1967` turns that
source into an explicit 2×2 matrix with a nonzero eigenvalue gap.  But an axis
is not yet a signed arrow: flipping the arrow still gives the same axis.  That
remaining sign is why global selector closure is still open.
