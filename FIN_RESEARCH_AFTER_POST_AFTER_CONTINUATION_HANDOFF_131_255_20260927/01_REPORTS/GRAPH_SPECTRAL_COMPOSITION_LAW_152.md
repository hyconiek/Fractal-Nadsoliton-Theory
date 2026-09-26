# GRAPH-SPECTRAL-COMPOSITION-LAW-152
## Once a graph coupling is declared, collective FIN gains are controlled by the graph Laplacian spectrum

Date: 2026-09-26

Status:
exact conditional theorem for PM-001.

Let

    m_x = X7^T p_x

and let L_Gamma be the weighted graph Laplacian of the declared cell graph.

The PM-001 quadratic retained interaction is

    (N/2)
    [
      g sum_x ||m_x||^2
      -kappa sum_(xy) w_xy ||m_x-m_y||^2
    ].

In matrix form:

    (N/2)
    m^T[
      (g I - kappa L_Gamma) tensor I_7
    ]m.

Diagonalize the graph Laplacian:

    L_Gamma v_r = ell_r v_r.

Then every graph normal mode v_r carries an effective retained FIN gain

    boxed:
    g_r = g-kappa ell_r.

## 1. Consequences

The uniform graph mode has

    ell_0=0,

so it always sees the original single-cell gain

    g_0=g.

Every nonuniform collective mode is shifted downward by

    kappa ell_r.

Thus the collective phase/bifurcation thresholds are predicted, CONDITIONALLY,
by

    g-kappa ell_r = g_threshold(single cell).

## 2. Examples

### two-cell edge

Graph-Laplacian eigenvalues:

    {0,2}.

Hence:

    g_common=g,
    g_relative=g-2kappa.

### three-cell path

Graph-Laplacian eigenvalues:

    {0,1,3}.

Hence collective gains:

    g,
    g-kappa,
    g-3kappa.

At the old scout point

    g=3.7,
    kappa=0.25,

these are

    3.7,
    3.45,
    2.95,

matching the previously reported path3 auxiliary minimum gain.

## 3. What this gives FIN

If a future law derives

    Gamma,
    w_xy,
    kappa,

then the collective retained-mode spectrum follows without additional fitting.

This provides a clean held-out test:
derive the graph/coupling first, then predict which relative modes soften.

## 4. What it does not give

The graph-spectral formula does NOT source the graph or kappa.

Different graphs with the same one-cell FIN core produce different collective
spectra.

PM-002 also shows that even the choice of retained-mean edge cost is not unique.

So this is a conditional propagation theorem, not a fundamental composition
law.

## 5. Next physically meaningful task

The next campaign should not fit kappa to a desired collective threshold.

It should search for an operational or deeper-relational quantity that:
- defines incidence;
- fixes edge weights;
- fixes the coupling normalization;
- survives a two-cell -> three-cell held-out test.

If no such source exists, FIN remains a well-controlled single-cell statistical
theory plus a conditional graph-coupling layer.
