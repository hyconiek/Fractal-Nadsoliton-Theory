# Second-pass audit of FIN-HYP-HIST-ME-001

Date: 2026-10-20. **Same-author second pass with a fresh computational
implementation; not an independently authored or external peer review.**
Original sources, model contract, self multiplicity, carrier and target are
unchanged. The full source proof and its controls were inspected rather
than taking the producer's PASS as a mathematical premise.

## 1. Audit dispositions

1. The stationary flow set is compact and nonempty. Positive diagonals
   allow a strictly positive feasible flow by assigning small symmetric
   off-diagonal masses and retaining the remainder on the diagonal.
2. Strict concavity gives one optimizing flow. The -t log t gain at a zero
   allowed entry, compared with order-t changes at positive entries, pays
   positivity. Transposition preserves both margins, so the optimum is
   symmetric. Reversibility was not assumed before this step.
3. The multiplier argument is valid despite one redundant marginal
   constraint: orthogonality to the constraint kernel suffices. Positive
   diagonals then determine z uniquely from c_xx=K_xx z_x^2.
4. The entropy-rate bound applies to **all stationary histories**. Equality
   forces all finite-past conditionals to be the one-step conditional and
   hence the same finite cylinders. No Markov assumption is smuggled in.
5. Simultaneous label translations and uniqueness make z constant on each
   difference fibre. Direct neighbour generation confirms one self versus
   two off-fibre continuations, including their multiplicities and 1/q
   Gibbs normalization. The strong quotient is valid, not a mean-rate fit.
6. The scalar-root equation and its lower/upper endpoints are correct.
   Nonuniform p amplifies next-destination odds, giving the claimed HB
   mismatch. The q>=3 and positive-p hypotheses are essential to that test.
7. The thinning double ratio, micro/macro entropy correction, marked-self
   non-adoption and separate binary N3 counterexample keep their stated
   scopes. Recovering HB after changing the reference is not a repair of
   the frozen primary hypothesis and does not supply a natural law.

No concrete mathematical counterexample to those scoped source statements
was found in this pass; no postulate or source file was repaired. This
disposition is conditional analytical admission within a self-review, not
closure of the promised independent-peer audit.

## 2. A different existence/uniqueness derivation: coercive dual

Let K be a fixed finite symmetric reference, positive on its allowed
support, with K_xx>0. Let pi_x>0 and sum pi=1. Set z_x=exp(s_x) and

    G(s)=1/2 sum_xy K_xy exp(s_x+s_y)-sum_x pi_x s_x.

Then

    partial_x G=z_x(Kz)_x-pi_x,
    v^T Hess G v=1/2 sum_xy K_xy z_x z_y(v_x+v_y)^2
                 >=2 sum_x K_xx z_x^2 v_x^2>0          (1)

for v!=0. The diagonal terms pay positive definiteness in every direction,
not merely along an assumed mass constraint. Moreover

    G(s)>=sum_x [K_xx exp(2s_x)/2-pi_x s_x].            (2)

Every summand is bounded below and diverges at both ends. In finite
dimension (2) therefore proves coercivity of G: no component can escape
to either infinity within a sublevel. G has a unique minimizer, its
gradient vanishes, and z_x(Kz)_x=pi_x. This proves the positive scaling
without using the source's entropy-boundary or multiplier construction.

The resulting c*_xy=K_xy z_x z_y has both margins pi and total mass one.
At the minimizer,

    G(s*)=1/2-sum pi log z,
    L_K*=2G(s*)-1+sum pi log pi
         =-2sum pi log z+sum pi log pi.                (3)

For F_K(c)=-sum c log(c/pi_x)+sum c log K, margins pi give

    L_K*-F_K(c)=KL(c||c*)>=0.                          (4)

Thus the dual gives the same optimum and exact gap as the primal source
proof. For the primary unmarked model K is 0/1, so the reference term
vanishes. Weighted K here only covers the already declared reference
extension; it is not a new choice of physical prior.

## 3. Hypothesis guards that must not be dropped

Without positive diagonal entries, arbitrary positive pi need not be
feasible: on the two-node off-diagonal-only graph the two margins require
pi_0=pi_1. Even at uniform pi its scaling z_0 z_1=1/2 is not unique.
The source contract includes self loops and therefore avoids both failures.

Entropy rate alone does not control initial mismatch for a nonstationary
preparation. For iid reference pi=(1/2,1/2), a process with deterministic
X_0=0 and iid future pi has the same asymptotic entropy rate as the
reference but initial TV=1/2. The new quantitative result consequently
retains exact stationarity and the common pi, not a free iid-parent transfer.
Even fixing every one-state marginal is insufficient without stationarity
of the history: take X_0 fair, X_1=X_0, and independent fair states from
time 2 onward. Its asymptotic entropy rate is log 2 and all state marginals
are fair, but its first pair has KL=log 2 against iid fair history.

No physical truth of maximum entropy, pi, support, intrinsic step or natural
reference follows from this mathematical audit. Strict F01 remains open.
