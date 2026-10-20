# ID-001 audit: native entropy, source budget and finite-dimensional corollary

Date: 2026-10-19. Conditional analytical mathematics, not formal verification,
independent peer review or a new dynamical/physical law. Originals unchanged.

## 1. The native model and clock are correctly bound

The actual operational intake supplies B=2(1p^T-I),
p_j=softmax_j[(g/2)A7_0j] and P_t=e^-2t I+(1-e^-2t)1p^T. Both labelled
copies attempt refresh at rate1. The factor2 is supplied model time, not
calibrated seconds. The actual state-change rate is
rbar=sum_i p_i 2(1-p_i), not the attempted refresh rate2. All fixed finite
g give p_i>0. Uniform p is used only in the explicitly g0 control.

This reuses the exact N2 quotient, not the older ST383 strict-heat matrix.
ST383's numerical entropy is in bits/symbol for its own supplied law; it
cannot be substituted into this N2/nats/common-time contract.

## 2. Fano joining inequality and the KS power identity

The two submitted block proofs are correct. Stationarity gives one error
probability at every coordinate; neither independent errors nor a Markov
readout is assumed. Conditional entropy chain rule bounds H(X^n|Y^n) by
nF_12(e), F_12(e)=h2(e)+e log11. The deterministic partition readout has
rate at most h_mu(T_delta). Full stationary-history joining is essential.

For completeness, the only time-scaling identity needed can be proved
directly. For any finite partition P and integer m, let
Q=join_(j=0)^(m-1) S^-j P. Its n-block under S^m is the mn-block of P
under S, so h(S^m,Q)=m h(S,P), giving h(S^m)>=m h(S). Conversely the
S^m n-block of any P is coarser than its S mn-block, so h(S^m,P)<=m h(S,P)
and h(S^m)<=m h(S). Thus h(S^m)=m h(S). Applied to T_(1/m), this pays
h_mu(T_(1/m))=H/m without an arbitrary-clock theorem or ergodicity premise.

The source finite-delta inequality, monotonic-branch inverse and its
liminf m e>=rbar consequence are therefore valid in the declared class.
They are standard entropy tools specialized to the native law, not a new
general Fano or KS theorem. A stronger consequence is separately derived
in SUBSAMPLING_ERROR_FLOOR.md and must not be attributed to the producer.

## 3. An explicit remainder payment

Write A=1-sum p_i², C=sum p_i(1-p_i)log p_i and t=1-e^-2delta. The exact
diagonal-inclusive identity in the source is correct:

    h_delta=-sum_i p_i(1-t(1-p_i))log(1-t(1-p_i))
             -tA log t-tC.

For delta<=1/4, |t-2delta|<=2delta²,
1-delta<=t/(2delta)<=1 and |log(t/(2delta))|<=4delta/3. The diagonal
Taylor error is at most t². Expanding with no diagonal deletion gives

    h_delta=rbar delta log(1/delta)+c(p)delta+R,
    c(p)=rbar(1-log2)-2C,
    |R|<=delta²[2A log(1/delta)+(14/3)A+2|C|+4].

This proves the source O_p(delta²|log delta|) with an explicit allowance.
Indeed A<=1 and |C|<=H(p)<=log12 give a uniform upper remainder bound over
the fixed12-simplex. This does not give a uniform positive rbar/error floor
as g becomes unbounded and p can approach a point mass. No p_min is invented.

For a full stationary clock-resolved CTMC path-space shift, exact sampled
readout and this divergence force h(T1)=infinity by the same power rule.
This is the known phenomenon discussed in the supplied primary PDF, not a
new FIN source theorem. Its regularized jump-sequence quantity in §3.2 is
explicitly a different definition. The relevant local PDF/text were read;
neither remote historical access nor original downloads are reconstructed.

## 4. Compact R26 means a conditional corollary, not the whole nadsoliton

For a fixed polynomial Hamiltonian vector field on R26 and an actually
supplied compact two-sided invariant K, L=sup_conv(K)||DF|| is finite.
Mean-value estimates and Gronwall give an e^L-Lipschitz time-one map on K.
An initial epsilon e^-L(n-1) net is an (n,epsilon) spanning set. Ambient
covering growth O(r^-26) therefore gives h_top<=26L and h_mu<=h_top for
each invariant probability on K. No compactness/invariant measure is derived
from unspecified signs; velocities are part of the full26-dimensional state.

This corollary covers the offered compact smooth source only if its clock,
measure, readout and stationary joining are really the proposed interface.
It does not establish finite entropy of the infinite-dimensional NS/RP PDE.
Conversely, infinite dimension or singularity is **not itself an escape**:
any source still satisfying finite H and the other measurable-flow premises
is covered. Noncompact/singular-limit or infinite-entropy alternatives need
their own paid objects, not a new label or a bath postulated to fit HB.

## 5. Preserve observation and physical boundaries

The source's same-p/different-clock countercontrol requires c>0 when
writing B'=cB as another irreducible CTMC. Its phrase “any c!=1” must not
be read as licensing negative rates. At a fixed delta=1/m, the suspension
countercontrol has h(T_delta)=h_delta and h(T1)=m h_delta by the same power
rule; this suffices for the sampled scope without an arbitrary-clock claim.

Full stationary joining is stronger than one-window TV, weak convergence,
two-time covariance or a spectral histogram. Positive record entropy is
not thermodynamic entropy production: this native chain satisfies detailed
balance and its stationary edge-ratio entropy production is zero.

The ontological existence of intrinsic change does not supply this external
comparison clock or make finite-H/finite-dimension premises properties of
the nadsoliton. The result excludes a precise effective source class. It
does not select a generator, physical time, geometry, units, QM or ToE.
