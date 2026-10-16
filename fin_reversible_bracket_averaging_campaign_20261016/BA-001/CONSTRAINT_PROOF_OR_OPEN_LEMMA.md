# BA-001 — flux proof and the remaining component lemma

**Proved exactly:** reversibility plus the statewise parity drift fixes the equilibrium traffic crossing every parity level, the parity-compressed instantaneous generator, and the whole parity-projected carré-du-champ. It does **not** fix traffic inside each parity fibre, the Markovian law of the projected path or the invariant components of the full microscopic chain.

For a single allowed edge `x<->y` put `c_xy=pi_x q_xy=pi_y q_yx`. All unknown activities are nonnegative and subject to the statewise constraints `sum_y c_xy(m_y-m_x)=pi_x d_native(x)`, the per-vertex cap `sum_y c_xy/pi_x <= CN` and reflected support. This is the precise linear finite feasibility problem before the open temporal question.

Let `A_l={x:m_x<=a_l}`. Sum the flux equation over `A_l`: each internal edge contributes `c_xy(m_y-m_x)+c_yx(m_x-m_y)=0`; same-parity edges also contribute zero. Every edge crossing the cut has jump `+delta` leaving `A_l`, so the sum equals `delta K_l`, as claimed. In particular, if some proposed deterministic diagnostic field sets `K_l<0` or differs from this identity it is not admissible.

One component `C` may be closed under the diagnostic rates even though many native edges cross its boundary; the same cancellation argument applies *within C* but involves its restricted Gibbs weights, not the unconditional parity marginal. Stationarity of `m²` on C gives `0=E_C[2m d + Gamma]`, hence the explicit formula for `h_C`. The logical gap is the uniform q12 family statement `sum_C pi(C) h_C² ->0`, or another additive-functional bound sufficient without this condition. Neither follows from the global parity-cut identities by algebra alone.

**Open, not asserted:** whether the specific q12/A7 MT geometry and exact native drift prohibit macroscopic component-by-component bracket variation even when some allowed activities vanish. Proving this requires a class-specific feasible-flux argument, not a proof for an arbitrary six-state graph.
