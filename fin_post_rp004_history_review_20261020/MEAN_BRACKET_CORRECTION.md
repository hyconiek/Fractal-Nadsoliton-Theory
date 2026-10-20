# Correction: a stationary mean bracket suffices in the reversible one-hot class

Date: 2026-10-20. New analytical correction derived in this intake. It is a
finite-class inverse identity, not an independent physical activity source.

## 1. The source sentence that is too broad

Section B of the incoming `INTRINSIC_HISTORY_PROVIDER_OR_OPEN_LEMMA.md`
correctly proves `q_ij=(Gamma_i)_{jj}` for the full conditional one-hot
increment tensor. Its following sentence says that an averaged quadratic
variation does not supply the conditional tensors. **That is false for
the full stationary mean matrix, known positive pi and reversibility in
this very one-label class.** A scalar average or only the diagonal is a
different, weaker observation. Do not silently import the many-count-state
AC insufficiency statement into the one-hot class.

The same overbroad statement is repeated in the source NONCONCLUSIONS.
Preserve the original files; admit their identifiability claims only with
the correction below. Their physical source-open conclusion remains valid.

## 2. Exact reconstruction and compatibility theorem

Let Q be a finite conservative jump generator, with known strictly positive
stationary probability pi and detailed balance `pi_i q_ij=pi_j q_ji`.
Observe the one-hot state `X=e_i`. Define

    Gamma_i=sum_(j!=i) q_ij (e_j-e_i)(e_j-e_i)^T,
    G=sum_i pi_i Gamma_i.

For distinct i,j, the only contributions to G_ij are i->j and j->i:

    G_ij=-pi_i q_ij-pi_j q_ji=-2pi_i q_ij.             (1)

Every increment has zero coordinate sum, so `G 1=0`. Its diagonal therefore
also gives the correct negative diagonal of Q. Writing `D_pi=diag(pi)`,

    G=-2 D_pi Q,
    Q=-(1/2)D_pi^(-1)G.                              (2)

Thus **the full stationary mean tensor alone uniquely recovers Q and then
every Gamma_i** in the declared class. No partition function or statewise
conditional measurement is needed once pi is independently known.

Conversely a real matrix G is such a stationary mean tensor for some
reversible generator with this pi iff

    G=G^T, G1=0, and G_ij<=0 for i!=j.                (3)

Necessity follows from the increment formula. For sufficiency set Q by (2).
The off-diagonal rates are nonnegative, row sums vanish, and detailed balance
holds since `D_pi Q=-G/2` is symmetric. Also `pi^T Q=0` by symmetry and
zero row sums. Formula (1) verifies that its stationary mean bracket is
exactly G. Strictly positive rates require strict off-diagonal negativity
on the declared complete support; a sparse support needs matching zeros.
Positive semidefiniteness follows from (3), but PSD alone is insufficient:
`G=u u^T`, `u=(1,1,-2)`, has zero row sums and is PSD, yet `G_01>0` gives a
negative proposed q_01 and is inadmissible.

For the uniform twelve-label radial comparison family, (2) reduces to

    G=-Q/6=A/6, q_ij=-6G_ij (i!=j).                  (4)

These equalities describe the chosen **one-label** generator, not the
interacting native N2 parent, occupation coordinates or a physical field.

## 3. Why the hypotheses and observation type matter

Without detailed balance the off-diagonal mean reveals only
`pi_i q_ij+pi_j q_ji`, not its antisymmetric current. For example on three
labels with uniform pi, take off-diagonal rates 1 in both directions, or
rates `1+kappa` along the directed 0->1->2->0 cycle and `1-kappa` in reverse,
with `0<kappa<1`. Every exit rate is 2 and pi is stationary in both cases.
Their full mean G is identical, while the oriented generators differ.
This supplies an exact countercontrol for dropping reversibility.

In the native many-copy count process, many microscopic edges at many
occupation states have the same coordinate displacement. A stationary mean
over those states sums multiple unknown activities into a single matrix
entry. The one-hot unique-edge argument (1) does not apply there. Existing
AC wheel/stationary-average obstructions remain unchanged.

If only the anonymous equality predicate is observed, the labelled entries
of G are not automatically available. Calling the matrix a relational
observable does not derive this observation algebra. Likewise units of a
per-time bracket require the declared comparison time. Neither is supplied
by an anonymous return function or a static equilibrium histogram.

## 4. Separate event-index corollary

For a stationary reversible discrete jump chain J, the mean **per-step**
one-hot increment tensor is

    G_step=2D_pi(I-J),
    J=I-(1/2)D_pi^(-1)G_step.                         (5)

This follows by the same two-edge calculation, retaining possible self
steps (which contribute zero). For the supplied radial true-jump class,
J_ii=0, so `G_step,ii=2pi_i`. When its constant CTMC exit rate s is supplied,
`Q=s(J-I)` and `G_time=s G_step`. Erasing holding times still leaves s
undetermined; (5) does not generate a Poisson law or physical seconds.

## 5. Consequence for the next source gate

An actually independent, labelled full **mean** bracket and pi would already
be a sufficient provider in this reversible one-hot class; demanding the
entire pointwise family is unnecessarily strong. This is a real relaxation
of a sufficient inverse-data interface, not closure of F01/F03. The present
package supplies G only as a function of already assumed Q, so it still has
no independent activity source. The conclusion remains conditional, and
cannot be transferred to a different carrier/readout without a new proof.
