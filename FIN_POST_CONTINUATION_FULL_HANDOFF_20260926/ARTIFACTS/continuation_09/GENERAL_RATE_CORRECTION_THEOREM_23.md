# GENERAL-RATE-CORRECTION-23 — first-moment robustness and drift tomography

Status: **PROOF_GRADE_WITH_DIRECT_OPERATOR_REPLAY** inside the declared stationary heat-bath + ME7 lane.

## Statement
Let an O(1/N) convention correction add a state-independent jump-rate matrix K=(K_ij) to the first correction of the finite-N empirical generator at the uniform state. Define its net first jump moment

a = sum_{i != j} K_ij (e_j-e_i).

For any real retained probe phi in V7 put

m = <phi^2>_u,   h = P_H(phi^2),

where P_H projects to the discarded k=1,2 sector. Then the change of the leading quartic ME7 closure coefficient is

Delta C_K(phi) = 36 m <P_H a, h>.

Thus the entire matrix K enters only through the four-dimensional hidden projection P_H a. Any correction with a=0 is invisible. Any correction with a in V7 is also invisible. This strictly contains the earlier target-matrix result.

## Operator proof
Write G0 for the uniform leading refresh operator, D for the full-minus-ME hidden-departure operator, and C_K for the O(1/N) convention operator. On a quartic F=(phi.d)^4 only compositions with one of each operator can contribute to the convention shift at the required order.

The degree-leading identities are

G0 (phi.d)^4 = 12 m (phi.d)^2,
D (phi.d)^2 = h.d,
D (phi.d)^4 = 6 (phi.d)^2 (h.d),
C_K f = a.grad f + terms lowering degree by at least two.

Because phi is retained and h is hidden, <phi,h>=0. Of the six permutations of G0,D,C_K, exactly three survive. Each equals 12 m (a.h); the other three vanish by retained/hidden orthogonality. Hence

Delta C_K = 3 * 12 m (a.h) = 36 m <P_H a,h>.

Terms of C_K below its first-moment derivative lower polynomial degree too far to survive the remaining two operators, so no higher jump moment of K enters this coefficient.

## Relation to the target-B theorem
For q^B_{j|i} proportional to exp[g(A7 p)_j-(g/N)B_ji], after column centering the induced rate correction has

a = -(g/12^2) B 1.

The general theorem therefore gives

Delta C_B = -(g/4) <phi^2>_u (B1)^T P_H(phi^2),

exactly the previously derived formula. In particular A7 1=0 proves exact leading-order invariance of Gibbs leave-one-out.

## Direct replay
A random completely nonsymmetric K was evaluated by the six exact finite-difference operator compositions. The three surviving permutations were equal to machine precision and the sum agreed with 36 m a.h to about 1e-14.

## Noncircular tomography of P_H a
Two translation orbits suffice. With orthonormal real Fourier modes define

phi^(1)=(c3+c4)/sqrt(2),       P_H[(phi^(1))^2]=(sqrt(6)/12)c1,
phi^(2)=c5,                    P_H[(phi^(2))^2]=(sqrt(6)/12)c2.

Translate each probe through all 12 labels. The intrinsic quartic coefficient is translation invariant and therefore occupies only DFT mode zero. The convention drift occupies DFT k=1 for the first orbit and k=2 for the second. Hence all four coordinates of P_H a are recovered from nonzero Fourier components without subtracting or assuming the intrinsic k=5 coefficient.

A random-K fixture recovered all four hidden coordinates with max absolute error 7.8e-16.

## Negative controls and scope
- Any K with zero net first jump moment is exactly invisible at this order.
- Any net drift lying wholly in the retained sector is exactly invisible.
- Convention-robust null probes satisfying P_H(phi^2)=0 remain null for arbitrary K.
- An arbitrary hidden net drift can shift, cancel, or reverse a non-null quartic fingerprint. Therefore the fingerprint is not convention-universal without either a balance premise or drift tomography.

This theorem concerns the declared finite-N heat-bath/ME7 generator. It does not source a physical clock, activity law, laboratory observable, scale, QW-2191, L_total, SM/GR or ToE closure.
