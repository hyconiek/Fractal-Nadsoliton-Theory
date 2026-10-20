# CF-002: from causal intensity defects to complete path-law error

Date: 2026-10-20. Self-contained finite-source theorem. Coupling, point-
process likelihood and entropy methods are reused, not claimed globally
new. This pays a quantitative source-interface gap left by exact ES closure.

## 1. Hypotheses and defects

Use the finite hidden source, readout and supplied preparation in
[CF-001](EXACT_CAUSAL_FILTER.md). Let P be its observed path law on [0,T],
with initial label law m. Let Q_path be the CTMC law with generator Q and
initial law nu. All source rates and target rates are finite. For the
entropy formulas assume nu_i>0 wherever m_i>0 and q_ij>0 for i!=j.

The exact visible intensities alpha are derived by the source filter.
On its current row i=D_(t-), define

    D_vis(T)=E_P integral_0^T sum_(j!=i)|alpha_ij(t)-q_ij| dt,
    C_vis(T)=E_P integral_0^T sum_(j!=i)Phi(alpha_ij(t),q_ij)dt,
    Phi(a,q)=a log(a/q)-a+q, Phi(0,q)=q.               (1)

Full-source upper defects replace alpha by a_ij(Z_(t-)); call them
D_src and C_src. Since Q's current row is observed, conditional Jensen
in CF-001 (4) yields `D_vis<=D_src`, `C_vis<=C_src`.
These are **whole-history** defect norms, not covariance mismatches or
tests against present-state functions only. Initial mismatch is not omitted.

## 2. Coupling bound, no positive-support assumption needed

The complete unmarked cadlag path laws satisfy

    TV(P,Q_path)<=min(1, TV(m,nu)+D_vis(T))
                <=min(1, TV(m,nu)+D_src(T)).          (2)

Construct a maximal initial coupling of the two labels. Preserve the source
conditional state law given its initial label. Until first label mismatch,
couple jumps to destination j with common intensity min(alpha_ij,q_ij),
source-only residual (alpha_ij-q_ij)_+ and target-only residual
(q_ij-alpha_ij)_+. The source filter follows its own deterministic equation
and updates only on its own observed jumps. After the first mismatch use
independent continuation. Bounded rates give a constructive dominating-
Poisson thinning construction on every finite horizon; zero-intensity
filter updates are never attempted. Each marginal retains its original law.

Until mismatch, its instantaneous hazard is sum_j|alpha_ij-q_ij|.
Compensation gives a mismatch probability bounded by its stopped expected
integral, then by the full P integral in (1). Add initial mismatch. Since
equality in this coupling is equality of complete paths, the coupling
inequality proves (2). No simulation and no independence of source
observations is assumed. A direct finite `(Z,target)` coupling similarly
proves the looser D_src bound while retaining hidden internal transitions.

## 3. Exact visible likelihood entropy

With the positive-support conditions above,

    D_KL(P || Q_path)=D_KL(m||nu)+C_vis(T).             (3)

Here is a direct derivation, not an unread change-of-measure theorem.
For a source history with states i0,...,in and ordered actual jump times,
the conditional survival between jumps is exp(-integral alpha_out) and
the jump density is alpha_ij. This follows either from the killed-filter
construction or its conditional hazard equation. Its path density divided
by the reference CTMC density has logarithm

    log(m_i0/nu_i0)
       +sum_jumps log(alpha_ij/q_ij)
       -integral_0^T sum_j(alpha_ij-q_ij)dt.          (4)

The reference supports all finite source histories; bounded rates ensure
no explosion. Under P, compensate the jump sum in (4) to get (3). For a=0
use a log(a/q)=0. Bounded a and positive finite q make these terms integrable
on a finite horizon. Truncating log away from zero and taking the limit
justifies compensation there as well, since a|log a| is bounded near zero.

Using `log x<=x-1` gives the directly computable upper bound

    C_vis<=C_src<=E integral_0^T sum_j (a_ij(Z_t)-q_ij)^2/q_ij dt. (5)

For completeness, the entropy-to-TV step follows by coarse-graining to any
event A (the log-sum inequality), then the binary entropy divergence has
second derivative at least 4 and minimum at equality. Thus
`D_KL(P||Q_path)>=2 TV(P,Q_path)^2`. Consequently

    TV(P,Q_path)<=sqrt([D_KL(m||nu)+C_vis(T)]/2),       (6)

and the same upper bound with C_src or (5) is sufficient. This is a
finite-horizon probability bound, not a physical validation or a guarantee
that a currently unavailable provider pays the RHS.

## 4. Derive the best constant Markov summary, without adopting it as a law

For a fixed observed source law define

    w_i=E integral_0^T 1_{D_t=i}dt,
    n_ij=E integral_0^T 1_{D_t=i}alpha_ij(t)dt.

Assume these w_i and off-diagonal n_ij are positive. Among reference CTMCs
with the same initial law and unrestricted positive constant off-diagonal
rates, the unique minimum path-KL generator is

    q*_ij=n_ij/w_i, q*_ii=-sum_(j!=i)q*_ij.           (7)

Indeed, subtract the cost for Q* from (3); all terms depending only on the
source cancel, leaving exactly

    D_KL(P||Q_path)-D_KL(P||Q*_path)
                        =sum_(i,j!=i) w_i Phi(q*_ij,q_ij)>=0. (8)

Equality requires every q_ij=q*_ij. This is an approximation/identification
theorem **given the source law**, not a maximum-caliber postulate sourcing
that law. For zero fluxes the corresponding boundary treatment is separate;
it is not needed for the positive ES example.

For a stationary source, w_i=T p_i and q* is its state-conditioned mean
transition rate. If its stationary pair fluxes are reversible, the mean
one-hot bracket satisfies `G=-2diag(p)Q*`, even if the process itself has
memory. Thus the earlier G inverse can identify the best Markov summary
without certifying full Markov closure. In ES, Q* is precisely the native
reset B for every epsilon. A different constant generator cannot remove
the positive residual history cost: it only increases path KL by (8).
No linear-in-T formula for the residual C_vis(T) is asserted.

## 5. Exact zero, preparation and horizon limits

When m=nu and q_ij>0, (3) vanishes iff alpha_ij=q_ij P dt-a.e. on active
rows, since Phi>=0 with equality only at a=q. This is an exact law-level
causal criterion. D_vis=0 yields the same closure by (2), even when target
support has zeros. Neither statement is licensed by pair-conditioned rates.

D_src=0 is sufficient, but D_src>0 is **not** a disproof of observable
closure: source knowledge can reveal pre-sampled future marks that the
observer does not know. CF-003 gives an exact example. C_src is likewise an
upper bound, not a necessary entropy payment of the observed process.

If a family is studied for growing T_N, the integrated defects themselves
must vanish on that growing interval (plus the corresponding initial term).
A per-time small error, a fixed-time limit, stationarity or matched moments
does not pay that requirement automatically. No finite-N approximation or
new resource theorem is inferred from these conditional inequalities.

For a native N2 target q_ij=2p_j, p_j>0 and all sums concern its twelve-state
difference law. A putative source for the 144-state parent, the count chain
or a field theory needs its own state/readout and defect bounds; these
carriers are not identified by notation.
