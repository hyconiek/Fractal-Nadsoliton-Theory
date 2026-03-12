# N466 Current First Strict Pair1 Diagonal Sector O(2)-Cut Criterion Theorem

Status: `N466_DISCHARGED_CURRENT_FIRST_STRICT_PAIR1_DIAGONAL_SECTOR_O2_CUT_CRITERION_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-12`

## Goal

After `N465`, one strict fact is now fixed:

```text
the frozen host A = K_total + m0^2 I is isotropic on pair1 and therefore cannot supply an O(2)-cut there.
```

So any strict `pair1` axis-selection / `O(2)`-cut attempt must come from **additional**
non-translation-invariant structure, most naturally from the **local diagonal sector**.

`N466` isolates the exact mathematical condition under which a diagonal sector can break `O(2)` on `pair1`.

## Strict-admissible inputs reused

1. `QW-2190/QW-2191`
   - the strict 12-site ring scaffold and the real Fourier pair `pair1 = span{c1,s1}`.
2. `N465`
   - host kernel + scalar floor cannot cut `O(2)` on `pair1`.

`N466` is otherwise a pure linear-algebra / trigonometry lemma.

## Setup (pair1 basis and a diagonal profile)

Let `n=12`. For `i ∈ {0,...,n-1}`, define:

```text
c1(i) := sqrt(2/n) * cos(2π i / n),
s1(i) := sqrt(2/n) * sin(2π i / n),
V1 := span{c1,s1} ⊂ R^n.
```

Let `D : R^n -> R^n` be a diagonal operator in the site basis:

```text
D = diag(d0, d1, ..., d_{n-1}).
```

Write the `pair1` restriction matrix in the basis `(c1,s1)`:

```text
M1(D) := [[a1, b1],
          [b1, d1]],

a1 := <c1, D c1>,
b1 := <c1, D s1>,
d1 := <s1, D s1>,
Δ1(D) := (a1 - d1, b1).
```

`D` breaks `O(2)` on `V1` iff `Δ1(D) ≠ (0,0)`.

## Theorem (closed form and exact criterion)

### Claim 1. Closed-form Δ1 in terms of the diagonal profile.

Using the identities:

```text
cos^2(t) - sin^2(t) = cos(2t),
2 cos(t) sin(t) = sin(2t),
```

one obtains:

```text
a1 - d1 = (2/n) Σ_{i=0}^{n-1} d_i * cos(4π i / n),
b1      = (1/n) Σ_{i=0}^{n-1} d_i * sin(4π i / n).
```

So:

```text
Δ1(D) = ( (2/n) Σ d_i cos(4π i/n) ,  (1/n) Σ d_i sin(4π i/n) ).
```
∎

### Claim 2. Diagonal O(2)-cut on pair1 is equivalent to a nonzero “mode-2” Fourier coefficient.

Define the complex mode-2 discrete Fourier coefficient of the diagonal profile:

```text
F2(d) := Σ_{i=0}^{n-1} d_i * exp( i * 4π i / n )  ∈ C.
```

Then, by Claim 1:

```text
Δ1(D) = ( (2/n) Re(F2(d)) , (1/n) Im(F2(d)) ).
```

Therefore:

```text
Δ1(D) = (0,0)   <=>   F2(d) = 0.
```

Equivalently, the diagonal restriction `D|_{V1}` is scalar iff the diagonal profile has no Fourier component
at frequency `2` (i.e. at phase step `4π/n`).

So a diagonal sector supplies a `pair1` `O(2)`-cut iff and only iff:

```text
Σ d_i cos(4π i/n) ≠ 0  or  Σ d_i sin(4π i/n) ≠ 0.
```
∎

## Consequence (what a strict “physical accelerator” would have to deliver)

1. A constant diagonal floor `m0^2 I` has `F2(d)=0` and therefore cannot break `O(2)` on `pair1`.
2. By `N465`, the strict frozen host `K_total + m0^2 I` is already isotropic on `pair1`.
3. Therefore any real `pair1` axis selection must come from an additional diagonal (or other local) sector
   whose profile carries a nonzero `mode-2` Fourier component.

This is a precise, checkable target: to get a strict-core `O(2)`-cut on `pair1`, the theory must export a
strictly sourced diagonal/local profile with `F2(d) ≠ 0` (or an equivalent non-translation-invariant
ingredient).

## What N466 does not prove

`N466` does not prove:

1. that the FIN canonical diagonal sector actually has `F2(d) ≠ 0`,
2. that any currently exported diagonal profile is strict-core derived,
3. strict-core theta export,
4. strict-core selector closure or `QW-2191` discharge,
5. ToE closure.

