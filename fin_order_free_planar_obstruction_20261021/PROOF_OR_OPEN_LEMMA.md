# Five ports exclude every circular order

Date: 2026-10-21. Conditional analytical theorem with an exact finite
certificate. Same-author proof, not independent review or formalization.

## 1. The imported necessary condition and sign convention

Use M=-L12, so M_ij=w_d>0 for distinct ports. Definitions 2.2.1/2.2.4 and
Theorem 2.2.6(a) in [Alman–Lian–Tran, arXiv:1309.2697v1](https://arxiv.org/pdf/1309.2697v1)
give the circular-minor condition in this sign convention. We use only its
necessity for positive disk networks, not a graph reconstruction theorem.
The relevant definitions, current-sign convention and theorem were read in
the actual PDF on printed pages 3–5. This does not claim an independent
check of the original Curtis–Ingerman–Morrow proof.

For four ports appearing cyclically as a,b,c,d, two permitted minors are

    det M[(a,b),(d,c)] = M_ad M_bc - M_ac M_bd >=0,
    det M[(b,c),(a,d)] = M_ab M_cd - M_ac M_bd >=0.           (1)

Thus the product for the crossing matching (ac)(bd) must be no larger than
either noncrossing matching product. If the three products have a unique
minimum, that matching is forced to be the crossing matching.

This statement concerns submatrices of the ORIGINAL twelve-port response.
We neither call a five-port principal submatrix a new zero-row-sum response
nor eliminate the other seven ports. Interspersing other ports does not
change the induced cyclic order of a chosen quartet.

## 2. Three forced matchings

The exact weight enclosures in CERTIFICATE.json prove each minimum below
is strictly smaller than BOTH alternatives, with every gap >1/20000:

| Quartet | Forced crossing matching | Minimum product | Other products |
|---|---|---|---|
| 0,1,2,3 | (02)(13) | w2^2 | w1^2, w1 w3 |
| 0,1,3,6 | (06)(13) | w6 w2 | w1 w3, w3 w5 |
| 0,2,3,6 | (03)(26) | w3 w4 | w2 w3, w6 w1 |

These statements are exact-real consequences of the frozen finite kernel,
not high-precision point comparisons or additional assumptions on ordering.

## 3. An order contradiction, without twelve-port enumeration

The first quartet fixes the induced order to 0,1,2,3, up to rotation and
reversal. Choose that orientation. In the second quartet, crossing (06)
with (13) forces port 6 onto the open arc from 1 to 3 containing port 2.
In the third quartet, crossing (03) with (26) forces port 6 onto the open
arc from 3 to 0, since port 2 lies on the other arc from 0 to 3.

Those open arcs are disjoint. No cyclic order of the five distinct ports
satisfies all three conditions. Therefore no cyclic order of all twelve
ports can make M the response of a finite positive disk network.

The finite checker also examines all 4!/2=12 five-port cyclic orders,
anchoring port 0 and quotienting reversal. It exports an actual negative
2x2 circular minor for each, with its ordered row/column lists. This is an
independent combinatorial check of the short arc argument, not a substitute
for the sign theorem. Every pair of the three constraints is satisfiable;
the displayed triple is irredundant as an order-constraint certificate.
No general minimality among all conceivable certificates is claimed.

## 4. End-to-end rational enclosure

For each d, put x=(743/4000)d+13/80. It lies in [0,13/10]. The cosine
Taylor terms decrease in absolute value, since x^2/2<1 and successive
factorial denominators increase. The degree-22 odd-index partial sum is
a lower bound; the degree-20 even-index partial sum is an upper bound.

For d^(9/5)=d*(d^4)^(1/5), use an exact integer fifth-root bracket at
denominator 10^24. Integer powers verify its endpoints. Positive interval
division encloses w_d. Outward rounding to twelve decimal places gives the
public weight boxes; their products and differences are computed as exact
Fractions. The six strict margins are verified from the WHOLE boxes.

No libm, numpy, mpmath.iv or unbound spectral estimate is used to decide a
sign. The previous natural-order certificate is preserved as prior evidence,
not silently promoted to this stronger conclusion. The model inputs remain
conditional: the rational calculation does not source these parameters.

## 5. Robustness to noncirculant perturbations

Let Mtilde be any symmetric candidate with each off-diagonal entry differing
from the frozen M by at most epsilon=1/100000. Its diagonal may be chosen
to make row sums zero; diagonal entries do not enter the displayed minors.
The original off-diagonal weights lie between zero and 1/2. A product of
two entries changes by at most epsilon+epsilon^2; a difference of two
products changes by at most 2epsilon+2epsilon^2.

Consequently each forced-match gap for Mtilde is at least

    delta=1/20000-2epsilon-2epsilon^2>0.                     (2)

The same contradictory crossing requirements remain. Mtilde need not be
circulant or preserve D12 symmetry. In every order at least one necessary
circular minor is <=-delta. Thus even this entrywise neighbourhood is
disjoint from the tested response class. Exact positive common scaling
also preserves the order-free obstruction; it is not a clock-unit issue.

In particular no sequence of such finite responses converges entrywise to
the supplied M, even if the port order changes along the sequence. Eventual
entrywise error <=epsilon would contradict the same excluded neighbourhood.
This excludes a proposed infinite/continuum realization ONLY IF its actual
boundary response is proved to be such a finite positive-disk limit. No
unproved discretization or general continuum-locality premise is imported.

## 6. What this establishes, and what it does not

The old failure depended explicitly on the natural port order. This theorem
removes that escape route for the supplied finite L12 and the stated
robust neighbourhood. It rules out an entire passive single-disk boundary
realization class, regardless of hidden graph size or chosen port order.

It does NOT exclude nonplanar networks, multiple boundary components,
distributed or dynamical media, vector/constrained extensions, other
operators, or all notions of locality. It does not construct a physical
bulk, determine spatial dimension, derive A7/kinetics or close F01. No
legacy/strict completion, selector or physical-role transfer is inferred.
