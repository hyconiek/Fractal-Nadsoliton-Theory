# IC: the clock-free part of the remaining F01 obligation

Date: 2026-10-20. New self-audited conditional model mathematics; not a
fundamental law, independent peer review or physical-source admission.
The ES source and CF filtering/likelihood theorem are reused, not rerun
as discoveries. The delta is a causal-clock quotient and an exact prepared
first-exit obstruction, including null refreshes. No RP task is executed.

## 1. Freeze the question before choosing a clock

Use an already supplied finite source/readout/preparation with bounded
observed-history intensities alpha_ij(t). Source existence is not derived
by this assumption. Fix a finite comparison generator Q with q_ij>0 for
i!=j and q_i=sum_(j!=i)q_ij. Put, on its actual current row i,

    a(t)=sum_(j!=i)alpha_ij(t),
    r_j(t)=alpha_ij(t)/a(t) when a(t)>0,
    s_ij=q_ij/q_i.                                      (1)

Here r is the conditional destination distribution **at an event**, given
the full observable past. It is not a stationary state-only average, nor
an unconditional transition matrix between consecutive observations.

Let Q^v denote a law with predictable intensities v(t,history)q_ij, the
same initial law, and a bounded nonnegative scalar speed v. It may depend
on the entire observed past and elapsed waiting time. It cannot choose
different speeds for different future destinations. A full causal source
filtration is also allowed when its q-branch martingale property is paid;
merely knowing a future mark and changing its clock is not that property.
Initial/source conditioning must not be silently changed.

This is a broad analysis comparison class, not a supplied physical clock.
An arbitrary anticipative reparametrization is distinguished in section 5.

## 2. Exact causal-clock entropy decomposition (IC-001)

For Phi(x,y)=x log(x/y)-x+y, extended in the usual nonnegative sense,
direct algebra gives

    sum_(j!=i) Phi(alpha_ij,v q_ij)
       =Phi(a,v q_i)+a KL(r||s_i).                      (2)

Set a KL=0 at a=0. If v=0<a both sides are infinite. When a>0 and v>0,
substitute alpha_j=a r_j, q_ij=q_i s_ij and sum r_j=s_ij=1;
the logarithm splits into log(a/(vq_i))+log(r_j/s_ij), proving (2).
The zero cases follow by limits. There is no entropy extremization premise.

CF's finite-source path-density calculation applies with the predictable
reference intensities v q_ij. Conditioning successive actual jumps gives
the same survival/jump product, now with integrated v. Consequently, with
the same initial law and the required support,

    KL(P_path||Q^v_path)
      = E_P integral_0^T Phi(a,vq_i)dt + C_branch(T),
    C_branch(T)=E_P integral_0^T a KL(r||s_i)dt.          (3)

The formula is extended-valued if a reference forbids a source event.
With bounded finite-source alpha, the minimizing predictable speed is

    v_*(t,history)=a(t,history)/q_i,
    inf_v KL(P_path||Q^v_path)=C_branch(T).              (4)

To pay attainment, compute alpha as CF's deterministic finite-source filter
functional on possible histories, and extend it boundedly to impossible
ones. Then v_* is bounded since min_i q_i>0. A common dominating Poisson
construction yields the comparison law. Wherever the source has a positive
jump intensity, a>0 and every reference destination has positive intensity.
The likelihood is thus supported. Its entropy integrand at the optimum is
bounded by a times a constant depending on min s_ij; a r_j log r_j has its
continuous zero convention. This pays the finite-horizon expectation, not
just pointwise optimization. A zero source total is matched by v_*=0.
If only strictly positive clocks are allowed, zero-total intervals give an
infimum by positive approximants; do not assert attainment there.

Hence the *history-dependent branching cost* cannot be removed by fitting
any such clock. It is invariant to a constant rescaling of Q. It is not
derived by calculating a stationary mean bracket. Conversely C_branch=0
only identifies the law with a causal-clocked Q in this prepared class;
it need not be homogeneous Markov in the original comparison time.

## 3. Exact closure in an operational clock, with its hypotheses

When r=s_i P dt-a.e. on a>0, the actual intensities are v_* q_ij. If the
accumulated clock A_t=integral_0^t v_*du is continuous, strictly increasing
and unbounded, its inverse supplies a Q-chain in operational time.
To prove this without importing an unverified clock theorem, condition on
the observed past just after an event. Until the next event, the state i
is fixed: integrated hazard is q_i times the increment of A. Thus the
next operational waiting time has survival exp(-q_i u), and the next
destination has probability q_ij/q_i independently of that waiting time
and earlier history. Iteration gives the finite Q path law. The assumptions
ensure the inverse exists through every finite operational horizon. Flat
clocks, finite total operational time or more general changes require their
own endpoint treatment; no such extension is silently used here.

Countercontrol: in two labels there is only one possible next destination,
so C_branch=0 for every finite intensity law. A renewal process alternating
labels with age hazard a(t)=1+t/(1+t) has survival (1+t)exp(-2t), not an
exponential. Its operational holding time uses A(t)=2t-log(1+t), with
survival exp(-A). Timing memory alone is not an irreducible branching
obstruction. This diagnostic is not a new FIN age variable or physical law.

## 4. Exact first **observed** exit in the already existing ES family

Use r=3, positive p, sum p=1, |h_i|<=1, sum p_i h_i=0 and |epsilon|<1.
The inherited block law is

    T(a,b,c)=p_a p_b p_c[1+epsilon h_a h_b h_c].         (5)

Independent blocks, uniformly mixed phase and the independent comparison
Poisson refresh clock lambda give ES. All two-time laws and present-state
increment brackets are already known to equal those of the reset target
Q_ij=lambda p_j, j!=i. The observed law is stationary/reversal symmetric;
the directed hidden source is not reversible. These are reused premises.

Start at a **stationary comparison-time origin**, conditional on D_0=i.
Let c_(n,j) be the probability that the next n refresh marks are i and the
following mark is j!=i. This is an analysis probability; null refreshes
are not observed data. Write x=p_i, H=1+epsilon h_i^3 and
B_j=1+epsilon h_i^2 h_j. Averaging the three starting phases yields

    c_(0,j)=p_j,
    c_(3l+1,j)=x^(3l+1) p_j H^l (2+B_j)/3,
    c_(3l+2,j)=x^(3l+2) p_j H^l (1+H+B_j)/3,
    c_(3l+3,j)=x^(3l+3) p_j H^l (2H+B_j)/3, l>=0.    (6)

Proof: proper block marginals are products p. A complete all-i block
contributes H; a complete final (i,i,j) block contributes B_j. For one,
two or three initial null refreshes, the three phase weights are
(2+B_j), (1+H+B_j), (2H+B_j). Each further three null refreshes insert a
complete all-i block and multiply by x^3 H. Formula (6) follows for all l.
It retains repeated/null outcomes rather than treating every refresh as a
state change. The complete-block probability x^3 H=T(i,i,i)<1, since all
other block words have positive probability. Therefore eventual observed
exit has probability one and the geometric series converges absolutely.

The prepared first-exit destination law is

    J_ES^first(i,j)=p_j + p_j/[3(1-x^3 H)] times
       {x(2+B_j)+x^2(1+H+B_j)+x^3(2H+B_j)}.           (7)

It is normalized over j!=i by the exit partition just proved. It is NOT
asserted to be a Markov jump-chain transition matrix for subsequent events:
starting at an event is a different preparation from starting at a stationary
time. No Palm/age conditioning has been discarded.

For the reset target J_0(i,j)=s_ij=p_j/(1-x). Equivalently,

    J_ES^first(i,j)/s_ij
       =1+theta_i[h_j+x h_i/(1-x)],
    theta_i=(1-x)epsilon h_i^2(x+x^2+x^3)/[3(1-x^3 H)]. (8)

This is an exact finite probability formula, not a small-time expansion.

## 5. A clock-independent event obstruction (IC-002)

Choose distinct i,j,k with h_i!=0 and epsilon(h_j-h_k)>0. On the frozen
initial preparation define the unmarked bounded statistic

    F_T=1_{D_0=i, first observed exit tau<=T}
                [s_ik 1_{D_tau=j}-s_ij 1_{D_tau=k}].     (9)

For **every causal-clock reference Q^v** of section 1, E F_T=0. Stop its
counting compensators at the first observed exit: the two contributions
are s_ik v q_ij and s_ij v q_ik and cancel pointwise. This proof does not
assume a constant/state-only clock, independent records or fixed exit time.

For ES, let G_(n+1)(T)=Pr(Poisson(lambda T)>=n+1), the CDF of the n+1-th
refresh time. The clock is independent of the analysis marks, so (6) gives

    E_ES F_T=p_i sum_(n>=0)G_(n+1)(T)
                              [s_ik c_(n,j)-s_ij c_(n,k)]. (10)

The n=0 term cancels. Every n>=1 term has the same positive sign, since

    s_ik c_(n,j)-s_ij c_(n,k)
      =epsilon h_i^2(h_j-h_k) x^n p_j p_k
          H^floor((n-1)/3) /[3(1-x)].                  (11)

In particular no unknown Taylor remainder is needed:

    E_ES F_T >= epsilon h_i^2(h_j-h_k) p_i^2 p_j p_k
              /[3(1-p_i)] * [1-(1+lambda T)e^(-lambda T)] >0. (12)

The range width of F_T is s_ij+s_ik, so layer-cake integration for a bounded
function, or its positive/negative parts, gives uniformly over all these v

    TV(P_ES,[0,T],Q^v_[0,T]) >= E_ES F_T/(s_ij+s_ik).    (13)

This concerns complete **unmarked** paths; tau and its destination are
measurable there. No attempt marks have been added to the observer.

There is also a purely order-based consequence: at T=infinity, (7)--(11)
give a nonzero first-exit destination contrast. Every order-preserving
onto time reparametrization starting at zero preserves that first observed
destination, even an anticipative one. Thus no such clock can turn this
prepared ES law into the reset law. This latter exact nonconjugacy does
not require a comparison unit. The **finite-T uniform inequality** above
does require the causal Q^v class: an anticipative clock could selectively
move destinations past T, so it is not covered by (13).

## 6. A rational control and the actual twelve-label scope (IC-003)

For p=(1/2,1/3,1/6), h=(1/2,-1/2,-1/2), epsilon=1/2, choose i=1,j=0,k=2.
Then H=15/16, s_ij=3/4, s_ik=1/4 and (7) gives exactly

    J_ES^first(1,0)=1883/2502, J_ES^first(1,2)=619/2502,
    J_ES^first(1,0)-3/4=13/5004.                       (14)

At lambda=2,T=1 the coefficient in (12) is 1/1728. Since the first five
terms of the series for exp(2) sum to 7 and all later terms are positive,
exp(2)>7. Hence every causal Q^v satisfies

    TV(P_ES,Q^v)>1/3024,
    inf_v KL(P_ES||Q^v)=C_branch(1)>1/4572288.          (15)

The entropy lower bound follows from Pinsker, proved in CF, not from a
false converse of Pinsker. All quantities here are small model fixtures,
not measured probabilities or a practical FIN N guarantee.

For the actual admitted N2 **difference quotient**, reuse
`p_d=softmax_d[(g/2)A7_(0d)]`, lambda=2. For positive nonzero g this p is
positive and nonuniform by the admitted nonconstant A7 row. The existing
ES choice h_d=p_d-sum_l p_l^2 is mean zero, bounded, nonzero and preserves
every permutation preserving p. For a nonzero, p-mean-zero h on at least
three positive-probability labels there exist i,j,k as required with
h_i!=0 and h_j!=h_k: take an
extreme i unless all remaining h agree, in which case take an element
of that repeated remaining class instead. The repeated class cannot be
zero by mean centering and positive p. If needed exchange j,k.
Thus for epsilon!=0 this same symmetry-compatible ES family has the
clock-independent obstruction, without choosing new kernel coefficients.

The formulas hold for arbitrary positive p, including twelve labels.
No full 144-state parent extension, interacting count-chain counterfamily,
critical limit, physical clock or new law is asserted. RS's reversible-source
exclusion remains intact; its assumption is not silently imposed on ES.

**Scoped non-identification corollary.** Let the supplied data consist of p,
all reset two-time kernels, the present-state increment tensors and the
inherited observed stationary/reversal/state symmetries. Within the tested
ES quotient class these data do not select even a history law **modulo
order-preserving onto clocks with a fixed origin**. Indeed epsilon=0 and
epsilon!=0 share those data, while (7)--(11) put their prepared first-exit
laws in different clock classes. A rule using only these common data cannot
be a valid unique source selector across that class. This is not a theorem
against every future formalization of FIN or a full-parent countermodel.

## 7. What this actually derives for F01, and where it stops

F01 is not merely a missing speed unit: for this existing signature class
it also needs an independently justified **causal branching/history rule**.
A new source cannot repair its destination odds by retuning any scalar
causal clock. One necessary target relation, when claiming a time-changed
Q readout, is alpha_ij q_ik-alpha_ik q_ij=0 on every active history; an
exact source-derived nonzero value instead pays a falsifier. This relation
is not inferred from a stationary covariance or stipulated as a new axiom.

Nothing in the current ontology, static kernel or supplied Gibbs law has
here selected that rule. The new entropy minimization is an analysis fit,
not a variational law of nature. Existing hypotheses do not license choosing
its minimizer as the nadsoliton's dynamics. F01/F03/F06 remain OPEN_SOURCE;
CSP is not activated. This bounded derivation ends at a clock-free missing
dependency, not another source inventory or an invented fundamental law.
Positive C_branch is also not thermodynamic entropy production or an arrow
of time: the ES observed history retains its inherited reversal symmetry.
