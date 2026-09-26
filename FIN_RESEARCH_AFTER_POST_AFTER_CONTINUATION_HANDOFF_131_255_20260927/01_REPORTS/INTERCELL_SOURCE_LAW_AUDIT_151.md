# INTERCELL-SOURCE-LAW-AUDIT-151
## Existing FIN relational/composition results do not source the inter-cell coupling

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, commit `ad15a909...`.

Status:
source-law audit against the existing PHY-002/003 composition campaign.

## 1. PHY-002 already identifies the missing primitives

The existing relational-unit result explicitly lists as NEW primitives required
for multicell physics:

- cell index set V;
- interaction support graph Gamma;
- edge weights w_xy;
- coupling kappa;
- relative feature-frame maps if frames differ;
- transition generator.

Therefore these objects were never derived from the frozen single-cell FIN
core.

## 2. PHY-003 is explicitly conditional

The declared retained-mean model PM-001 is

    P
      proportional to
    prod_x multinomial(n_x) 12^(-N)
    exp[
      (Ng/2) sum_x ||m_x||^2
      -(N kappa/2)
       sum_(xy in E)
       w_xy ||m_x-m_y||^2
    ].

The report itself states that

    Gamma, w_xy, kappa

and the common feature frame are modelling premises, not FIN theorems.

It also constructs PM-002, a full-p coupling, with the SAME isolated cell and
zero-coupling limit but different coupled predictions.

For the two-cell N=2 scout the two admissible models differ by total variation

    0.04540446.

So composition existence was proved, while composition uniqueness was already
refuted.

## 3. Connection to report 149

For a two-cell PM-001 edge with unit weight, write

    m1=m+d,
    m2=m-d.

Then the exponent's quadratic interaction becomes

    g ||m||^2
      +(g-2 kappa)||d||^2.

Thus:

    common retained mode sees
      g_common=g;

    relative retained mode sees
      boxed:
      g_relative=g-2 kappa.

This is exactly the same one-parameter freedom found abstractly in report 149.

The free ratio can be written either as

    eta

or, in PM-001 normalization,

    kappa/g.

## 4. Sensitivity at the current working gain

At

    g=5.145228719489142

the k6 threshold is

    g6=5.123427551398616.

The two-cell relative mode can cross that threshold only if

    g-2 kappa > g6,

i.e.

    boxed:
    kappa < 0.010900584045.

So an inter-cell coupling as small as about

    0.2119% of g

changes the qualitative relative-phase content at this working point.

The previously used PM-001 scout value

    kappa=0.25

is far outside that narrow relative-k6 window.

This is not a problem for the scout; it demonstrates that the collective
physics depends decisively on a parameter the single-cell theory does not
source.

## 5. Source-law verdict

No target-blind inter-cell source law is currently present in the accepted
single-cell/refinement/relational results.

The refinement-conductance work supplies a COMPOSITION ALGEBRA once an edge is
given.

It does not determine:
- which cells are incident;
- which graph is physical;
- the edge weight;
- the coupling scale.

Therefore P0-151 returns:

    boxed:
    NO UNIQUE INTERCELL SOURCE FOUND.

This is a real negative result, not a search failure:
the older PHY campaign itself explicitly records the necessary quantities as
new primitives and exhibits incompatible admissible couplings.
