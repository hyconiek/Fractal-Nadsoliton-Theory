# D3-ORBIT-SPLITTING-67
## Four D3 parent states generate a twelve-state daughter orbit

Date: 2026-09-26

Status:
- exact group-orbit counting given the MP7-039 D3 isotropy and the nonlinear
  classification of report 66;
- local bifurcation statement only.

## 1. Parent isotropy

The two-harmonic MP7-039 base state is fixed by

    <T^4,R> ~= D3,

of order 6.

The full internal symmetry group D12 has order 24.

Therefore the full D12 orbit of the parent has size

    |D12|/|D3|
      = 24/6
      = 4.

So there are four symmetry-equivalent parent states in the full label orbit.

## 2. Local daughters

The standard D3 critical representation has three reflection axes.

Because the cubic coefficient is nonzero, report 66 finds three local daughter
rays on either side of the crossing.

A generic daughter on one reflection axis retains only that reflection:

    stabilizer ~= Z2,

of order 2.

Hence its full D12 orbit has size

    |D12|/2
      = 12.

## 3. Consistency of the counting

Each of the four parent states has three local daughter directions:

    4 x 3 = 12.

This is exactly the size of one D12 orbit with reflection stabilizer.

So the local D3 bifurcation naturally reorganizes

    boxed:
    4 parent states
      -> 12 daughter states.

On the two sides of g*, the two triplets are rotated by pi/3 in the local
critical plane, giving two locally distinct twelve-state daughter-orbit
continuations.

## 4. Why this is useful

This gives a strict group-theoretic prediction for numerical continuation:
near the crossing, a complete daughter census should appear in multiples
consistent with a 12-member D12 orbit.

## 5. Emergence interpretation

The mechanism is:

    higher stabilizer D3
      -> critical two-dimensional representation
      -> cubic invariant Re(z^3)
      -> three local choices
      -> lower stabilizer Z2
      -> twelve-state full orbit.

This is a clean mathematical example of

    symmetry -> instability -> discrete multiplicity.

It is potentially relevant to FIN's emergence programme, but no identification
with observed physical multiplicities is licensed.
