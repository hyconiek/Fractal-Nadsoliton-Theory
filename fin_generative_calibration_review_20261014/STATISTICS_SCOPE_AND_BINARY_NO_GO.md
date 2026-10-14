# Conditional statistical test and the fixed-binary-event obstruction

Date: 2026-10-14. Scope: the declared N2 static FIN/flat law classes;
zero physical observations. Exact probability statements below are conditional
on their sampling/model premises. Input certificates assume mpmath.iv.

## 1. What the finite-data theorem establishes

For iid observations of a thirteen-symbol alphabet, including INVALID,
let q be their common law and qhat their empirical histogram. For probability
vectors, TV(qhat,q)=max_S[qhat(S)-q(S)]. Each proper nonempty event has the
one-sided Bernoulli Hoeffding bound exp(-2n eps²). There are8190 such events,
so

    P(TV(qhat,q)>eps) <= 8190 exp(-2n eps²).

For an independently frozen hypothesis set H, reject only when a proved
lower bound on inf_(p in H)TV(qhat,p) exceeds eps. If q is in H and the
concentration event holds, H cannot be rejected. This gives conditional
false rejection <=alpha_validation. It does not prove calibration coverage.

If two sets have separation rho>2eps, an exact-distance test excludes the
wrong alternative on this event. A conservative lower-distance program may
still return inconclusive. “FIN_ONLY” in the prototype means **FIN not
rejected and flat rejected**, not proven membership of FIN or proof that
the finite family describes nature. The true law may fit neither model.

No scalar empirical ESS establishes the iid premise. Nor do reset IDs by
themselves establish independent initial conditions, stationary output,
constant detector parameters or absence of common drift.

## 2. Sample count and additional errors

At alpha_validation=.025 and eps=rho/3, independent interval evaluations
verify the following integer thresholds for this particular union bound:

| Conditional remaining rho | Sufficient independent records |
|---|---:|
| .00119114 | 40,278,676 |
| .00039114 | 373,539,331 |

For each row,8190 exp(-2n eps²)<.025 whereas the n-1 expression is greater
than .025. See verification.json. These are thresholds for a conservative
bound, not necessary sample sizes for every possible test or achievable
apparatus rates.

The second rho is the source F03's illustrative additional allocation:
.0002 burn-in+.0001 field/compilation+.0001 drift **per hypothesis**.
It consumes .0008 of the pair gap. The40-million design must not be combined
with this allocation while still claiming its original separation/power.
Revise the class, radius and sample plan before validation. None of these
numbers is measured calibration.

The total .05 error budget additionally needs a real simultaneous calibration
coverage theorem at alpha_calibration=.025. Merely reserving that probability
does not pay it. Properly disjoint calibration/validation data and a frozen
law inclusion event are required.

For fixed positive heat-bath p, worst-case initial TV is1-min(p), so the
exact worst-initial-law burn-in threshold is
`max(0, log((1-min(p))/eps_burn)/2)` in the supplied unit clock. The source's
`log(1/eps_burn)/2` is a useful conservative sufficient payment, not a
necessary waiting time for that fixed distribution or a physical clock.

## 3. Gain-mesh exclusion: theorem versus program

For p_g=softmax(g a/2),

    ||d p_g/dg||_TV = (1/4) E_(p_g)|a-E_(p_g)a|
                   <= range(a)/8.

The mean absolute deviation of a variable in an interval is at most half
its range. For a grid of width h, a nearest node lies within h/2. Therefore
the grid minimum distance minus range(a)*h/16 is a valid lower bound on
the continuous-gain minimum. The source uses twice that mesh payment,
range(a)*h/8, which is conservative in exact arithmetic.

For a probability-simplex TV ball of radius r around a nominal p,

    inf_(TV(v,p)<=r) TV(qhat,v) = max(0,TV(qhat,p)-r).

The lower bound is the triangle inequality; equality follows by moving from
p toward qhat along their convex segment. Hence subtracting the radius in
the source's exclusion rule is valid for its broad **output-law** class.

The binary64 grid minimum and subtraction are not outward evaluated, so
the source is not a certified distance oracle on arbitrary threshold inputs.
The generous mesh factor is not a quantified floating error bound. Before
deployment, reject only using validated arithmetic or a proved rounding
payment. Validate finite integer n, finite normalized counts and alpha.

## 4. All single fixed binary events fail uniformly in this class

Consider the two source classes

    H_X = union_(g in [2.97,3.03]) {q: TV(q,p_X(g))<=.005},

with separately allowed gains and output laws, thirteen outputs and zero
nominal invalid mass. Full-law separation remains at least .00119114 by the
admitted gain certificate. Nevertheless:

> For every one fixed event S, there exist q_FIN in H_FIN and q_flat in
> H_flat with q_FIN(S)=q_flat(S).

Thus no test observing only repeated indicators of that event can uniformly
distinguish the classes, regardless of sample count.

The intake checks every event. For each model it uses certified nominal
probability boxes containing the point gains2.97,3 and3.03. For each subset
of the twelve valid labels, choose one admissible gain pair. The interval
upper bound on the event-probability difference is **less than .00818** for
every selected pair, hence strictly less than the two .005 radii combined.

For each probability vector p the projection of its TV ball onto event
probability is [max(0,p(S)-r),min(1,p(S)+r)]. Every value in that interval is
attainable by transferring mass between S and its complement. The two
Bernoulli intervals therefore overlap, proving the statement.

The serialized certificate covers all4096 valid-label masks. Including or
excluding INVALID leaves the nominal event probabilities unchanged; the
empty/full valid masks supply the two boundary events. Consequently all
**8190 proper nonempty thirteen-output events** are covered. The exact-
fraction checker recomputes each selected upper bound and rejects omitted
coverage. See BINARY_EVENT_NO_GO_CERTIFICATE.json and
check_binary_certificate.py; no new transcendental producer is trusted here
beyond the already audited E probability enclosures.

The quantifiers matter: a **different** admissible law pair may witness each
event. This does not assert one pair agrees on all events or contradict
full-law separation. It strengthens the source's merely inconclusive
one-event estimate for the declared class.

## 5. Limits and permitted escapes

- This is a no-go for one fixed binary event, not for the full histogram,
  multivalued/weighted statistics, multiple recorded events, sequential
  protocols or an independently calibrated narrower nuisance class.
- It uses output-TV balls. Do not transfer it to an imposed shared known
  channel, a tighter per-state readout constraint or known gain without
  rechecking the class.
- Calibration can change the class and permit a cheaper event test. Such
  a test needs a new prospective proof; validation-based event selection is
  not licensed by this result.
- Positive identifiability and a conditional confidence theorem do not
  provide a physical apparatus, clock, independent A7 source or evidence
  that nature realizes any of the specified models.
