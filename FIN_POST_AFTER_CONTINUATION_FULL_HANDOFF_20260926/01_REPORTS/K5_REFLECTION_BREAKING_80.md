# K5-REFLECTION-BREAKING-80
## A subcritical reflection-breaking event connects the k=5 reflection branch to the generic g=5 atlas component

Date: 2026-09-26

Status:
- stationary locations and branch matching are NUMERICAL;
- symmetry and local index bookkeeping are exact once the numerical simple
  zero is accepted;
- no interval certificate yet.

## 1. Reflection parent branch

Consider the reflection-symmetric k=5-related component containing R7P-037
orbits

    i=4, index 1
and
    i=13, index 2,

joined by the fold

    g≈4.890011279965.

A representative reflection has a three-dimensional fixed subspace and a
four-dimensional odd complement.

## 2. First odd-sector crossing

On the i=13 side, an odd reflection-breaking eigenvalue vanishes at

    boxed:
    g_RB1 ≈ 5.152672504944.

At the crossing the full H7 spectrum has the local sign pattern

    two negative,
    one zero,
    four positive.

The zero eigenvector is reflection-odd.

Below the crossing the reflection parent has index 2; above it the parent has
index 3.

## 3. Generic daughter

Perturbing along the odd zero mode on the daughter side and performing
pseudo-arclength continuation produces a branch with trivial stabilizer.

Following that branch downward in g gives:

    generic daughter
      -> at g=5 matches atlas orbit i=14
         to numerical D12 distance < 1e-12
         with index 3
      -> fold at
         g≈4.913081079101
      -> returns through g=5 as atlas orbit i=5
         with index 2
      -> continues to larger g.

Thus component C from report 76 is not independent of the k=5 reflection
family.

It is born in a reflection-breaking event from component B.

## 4. Local type

The parent symmetry at the crossing contains a Z2 reflection.

The critical mode is odd under that Z2, so the reduced potential is even in
the symmetry-breaking coordinate y:

    Phi_red
      =
      Phi0
      +(lambda/2)y^2
      +(beta/4)y^4
      +....

The generic daughter exists on the side where the parent still has the smaller
Morse index.

Together with the observed indices

    parent:   2
    daughter: 3

this identifies the local event as a SUBCRITICAL reflection-breaking
pitchfork-type bifurcation.

The terminology refers to the local Z2 quotient, not to the earlier D3 event.

## 5. A second tiny angular crossing

The small-amplitude reflection sheet close to the uniform k=5 threshold has
another numerical odd zero near

    g_RB2 ≈ 5.228426082016.

Its amplitude is much smaller and its angular stiffness is of the expected
O(r^10) scale.

This second crossing is above the uniform k=5 threshold and belongs to the
delicate degree-12 angular regime of report 79.

It is not needed to explain atlas i=5/i=14 at g=5, and its daughter component
has not yet been globally assigned.

## 6. Reflection-sheet fold above the uniform threshold

The same reflection component also possesses another numerical fold near

    g ≈ 5.235454792465.

At that fold the reflection-fixed branch turns in g while several angular
eigenvalues are already extremely small.

This gives the near-uniform branch cluster a three-event structure:

    uniform k5 threshold
      g≈5.22055480

    secondary odd crossing
      g≈5.22842608

    reflection-sheet fold
      g≈5.23545479.

The closeness of these events is a consequence of the very weak degree-12
angular anisotropy.

## 7. Branch-graph consequence

The known g=5 components now satisfy

    reflection component B
      |
      | reflection-breaking at g≈5.15267250
      v
    generic component C, index 3
      |
      | through atlas i=14 at g=5
      v
    generic fold g≈4.91308108
      |
      | index 2
      v
    atlas i=5 at g=5
      -> high-g continuation.

So two pairs that looked unrelated in the raw atlas are parts of one larger
symmetry-breaking network.

## 8. Next atom

`K5-ANGULAR-EVENT-CERTIFICATION-81`

The most valuable certification target is the first reflection-breaking point
g_RB1, because it directly links two g=5 atlas components.

Use:
- the exact reflection fixed/odd decomposition;
- three parent stationary equations;
- one odd-block zero-eigenvalue equation;
- interval transversality and quartic-sign test.

Acceptance:
- interval-certified simple Z2 crossing;
- sign of the reduced quartic coefficient;
- rigorous local connection of the reflection and generic branches.
