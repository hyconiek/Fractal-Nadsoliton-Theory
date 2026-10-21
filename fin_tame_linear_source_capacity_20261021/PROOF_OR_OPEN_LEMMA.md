# LS: continuous seeds do not give a linear canonical readout unlimited histories

Date: 2026-10-21. Conditional mathematical source-class result, not a
physical no-go or an independent review of the imported theorems.

## LS-001. The published mathematical inputs and their exact roles

The real field with its usual exponential is o-minimal: this is Wilkie's
theorem, read in its author's [Theorem 5.2, printed p140](https://www.numdam.org/item/AST_2009__326__131_0.pdf#page=11).
The 2009 author exposition reports the original1996 result; the full
original proof was not independently rechecked here. No assumption that
the physical universe is o-minimal is made.

[Aschenbrenner–Dolich–Haskell–Macpherson–Starchenko, Theorem1.1](https://arxiv.org/pdf/1109.5438v1#page=4)
states the uniform shatter bound pi_phi(B)=O(B^n) for a fixed formula in
a weakly o-minimal theory with n parameter variables. Their definitions
and duality, printed pp2–4 and16, distinguish object and parameter roles.
We use this mathematical theorem, not claim a new VC theorem or that the
small verifier proves either imported result.

For our purpose write the formula as

    psi(A;xi): sum_(l even) exp(a_l dot xi+b_l)
                 >sum_(l odd) exp(a_l dot xi+b_l),          (1)

where A contains the twelve affine coefficient arrays. Crucially A are
the paper's OBJECT variables and xi are its PARAMETER variables, d of
them. For fixed queries A_1,...,A_B, varying xi cuts out subsets of this
finite set of query points. These subsets are exactly the source's output
bit words. This is the dual orientation relative to treating logits as
parameters of a classifier on initial-state space.

The formula uses only field operations, order and real exp. Hence the
two imported results give a finite constant C, depending only on this
fixed formula/q/d, such that

    number of output words on B fixed queries <=C(B+1)^d. (2)

Repeated queries cannot increase the count. C can be enlarged to cover
all finite B>=1. It is uniform over every real coefficient value, every
time choice and N-dependent parameter array; no coefficient magnitude,
initial-density regularity, compactness or bounded precision is assumed.
Its value is not computed, so no usable finite-N threshold is exported.

Strict > in (1) makes ties zero. Definability/counting includes all real
points, not only full-dimensional open regions. It therefore permits
atomic probability on ties. Arbitrary fractal/nondefinable initial support
or probability does not change (2): the set of all possible words was
bounded before choosing any initial measure.

## LS-002. The actual canonical map and accessible entropy

The repository's dual7 in fin_rank7_followup/src/model.py computes
h=X7 theta and p=softmax(h). With sigma_l=(-1)^l,

    1{sigma dot p>0}
      =1{sum_even exp(h_l)>sum_odd exp(h_l)}.              (3)

The normalizing denominator is strictly positive. If a proposed source
propagates theta affinely from a fixed-dimensional initial vector xi,
its bit is precisely (1). This is a comparison of a source to the actual
state-map type, not an adopted projection from a physical field to counts.
An arbitrary common logit shift cancels from (3) and does not supply extra
word capacity at the prescribed sample times.

Allow an additional initial discrete variable S_N of cardinality R_N,
with branch-specific coefficient arrays and any JOINT law of (S_N,xi).
The constant in (2) is uniform over branches. Union of their word sets gives

    number of possible source words <=R_N C(B+1)^d,
    H(Yhat_0,...,Yhat_(B-1))<=log R_N+log C+d log(B+1).    (4)

This is a bound on ACCESSIBLE RECORD entropy, not on the total entropy of
the continuous initial state. The latter need not have any finite entropy
or finite precision. Correlated label/field initial data are allowed.

For example,12^N initially labelled states plus a fixed number of real
field coordinates still satisfy (4). A static continuous auxiliary field
is therefore not a source of arbitrarily many random updates merely
because its initial value has infinite precision. Proper nonlinear or
noisy feedback dynamics is not excluded by that observation or by this
restricted source class.

## LS-003. Actual native finite-record consequence

Use BI's actual iid-parent bits at model times tau_N+sqrt(N)jt, fixed
critical t>0, with M_N=N^2 transitions and B_N=M_N+1 bits. Its paid native
entropy/coupling theorem gives

    H(Y_0,...,Y_M)>=M c_N(t)-epsilon_N B log2-b2(epsilon_N),
    c_N(t)->c(t)>0, epsilon_N->0.                         (5)

The horizon sqrt(N)N^2 t is polynomial, inside the admitted TR window.
No target trajectory is postselected. For ANY joint law with a proposed
source record, write e_N for its mean bit mismatch. Conditional chain rule
and binary Fano, without independent-error or stationarity assumptions, give

    H(Y_0,...,Y_M)<=H(Yhat_0,...,Yhat_M)+B b2(e_N).        (6)

Substitute (4) into (5)–(6):

    (B/M)b2(e_N)>=c_N(t)
      -[log R_N+log C+d log(B+1)]/M
      -epsilon_N(B/M)log2-b2(epsilon_N)/M.                (7)

If R_N<=12^N and d is fixed, the capacity term divided by N^2 vanishes.
Consequently

    liminf_N b2(e_N)>=c(t),
    liminf_N e_N>=b2^(-1)(c(t))>0, inverse on[0,1/2].     (8)

The lower branch is valid even if some actual mismatch exceeds1/2; a
subsequence below the stated floor would contradict continuity of b2.
This rules out vanishing mean-bit error for the complete source class,
even with arbitrary continuous initial probability and arbitrary parameter
fitting. It does not select a fundamental law or prove indeterminism.

No finite-H KS bound, invariant source measure, full stationary joining,
computability, stable source spectrum or physical clock calibration was
used. It is an information-theoretic statement about the specified model
experiment, not an efficient laboratory distinguishing protocol.

### Fixed discrete carrier: even weak matching of all fixed blocks fails

If R_N is bounded uniformly, source entropy is <=log R+log C+d log(M+2).
For a sufficiently large but FIXED M this is smaller than M c(t).
If its law were TV-close to the actual target at every fixed block length
as N grows, entropy continuity on the fixed finite alphabet would contradict
(5). Thus this fixed-size subclass cannot recover all fixed finite bit
laws either. This is not extended to growing discrete carriers such as
FH; those may match fixed windows while paying growing initial resources.

## LS-004. Linear amplitudes and waves: a scoped corollary, not Born physics

For fixed real dimension and fixed polynomial degree, a polynomial
threshold is another fixed formula in the ordered real field. Classical
quadratic-intensity readout is one such case. Applying the SAME theorem
with its initial coordinates as
parameters again gives O(B^d) words. This covers finite-dimensional
unitary/heat propagation and finite-mode linear waves with such a readout,
provided all real/imaginary and position/velocity coordinates are retained.
It does not treat Born-sampled outcomes as deterministic thresholds.

Time dependence may be oscillatory: at predetermined times, propagation
matrices are arbitrary fixed real parameters in each query. Their time
dependence need not itself be definable. A state-dependent internal clock,
random operator or feedback must instead be tested in its actual composed
source-to-bit map; it is not silently replaced by fixed query arrays.

A finite-dimensional random preparation of a LINEAR PDE can likewise be
covered if its proposed probes/logits depend affinely on those same finite
initial coefficients. This does not exclude an infinite-random-mode field.
A fixed-order polynomial perturbative truncation may be covered as a fixed
formula; its bit predictions do NOT automatically approximate an exact
nonlinear PDE merely because the field-norm remainder is small. A paid
threshold-margin/bit-error transfer would be required.

## 5. Why the hypotheses cannot be dropped

One real coordinate can encode arbitrarily many random bits. Nonadaptive
queries sin(2^j pi x)>0, or the nonlinear doubling map with a fixed binary
readout, access such digits. Their whole history family is not (1) or a
fixed tame initial-state formula with freely chosen query coefficients.
Definability of EACH step of an iterated map does not imply uniform
definability of ALL iterates. Adaptive binary search is another exception:
its next threshold depends on earlier answers, so query points are not
fixed before varying xi. The verifier exhibits both exceptions, not new
FIN mechanisms or observers.

The actual QW/legacy action has independent nonlinear self/cross terms.
It is not the linear class above; reduction to a single radial oscillator
would require an explicit coefficient theorem. No source integrability, chaos, physical role or
nonlinear-history capacity is inferred from a useful special coefficient
choice. Any proposed escape must supply its real composed map and law.

Growing source dimension, growing readout/formula complexity, nonlinear
history feedback, additional stochastic outcomes or an unbounded initial
discrete resource can leave the class. They are NOT proved adequate by
being named. Exponential accessible word capacity is necessary for the
native records, not sufficient to select their probabilities or source HB.
F01, ontology, physical units/clock/locality/conservation/calibration/EV,
selector/legacy-role and QM/GR/SM/L_total/ToE gates remain unchanged.
