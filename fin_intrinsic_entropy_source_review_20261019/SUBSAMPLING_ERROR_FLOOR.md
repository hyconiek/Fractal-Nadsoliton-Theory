# Additional intake theorem: a strictly positive stationary-joining error floor

Date: 2026-10-19. This is an additional analytical consequence of the exact
source contract, not a result retrospectively inserted in the submitted
ID-001 package. No new source, clock, observation law or premise is added.

## 1. Decimation retains the same error

Fix delta=1/m and any of the source's allowed full stationary joinings.
For every integer k>=1, take coordinates (X_(kn),Y_(kn)). This is a stationary
joining under the k-step shift. Its mismatch at coordinate0 is still e_delta.
Its target is exactly the N2 skeleton at tau=k/m by the semigroup law;
its source is the SAME f_delta observed under T_(k/m). The power rule gives
h_mu(T_(k/m))=kH/m. No consistent joining across different m is required.

Applying the already paid block argument to this decimated pair yields

    F_12(e_(1/m)) >= h_(k/m)-(k/m)H, for EVERY integer k>=1.       (1)

Thus a stronger finite-grid necessary bound is

    e_(1/m) >= F_12^-1(D_m(H)),
    D_m(H)=sup_(k>=1) [h_(k/m)-(k/m)H]_+ .                       (2)

The inverse uses the original branch0<=e<=11/12. This is a necessary
constraint, not an achievable rate-distortion function or optimal source.

## 2. Nonvanishing error as the grid refines

Define the target-dependent number

    D_*(H)=sup_(tau>0)[h_tau-tau H]_+ .

It lies in[0,log12]. For every finite H it is STRICTLY positive, since
h_tau/tau diverges as tau decreases to0 and rbar>0. Fix any tau>0 and
choose integers k_m=max(1,floor(m tau)); k_m/m tends to tau. Exact native
entropy is continuous in positive tau. Equation(1) therefore implies

    liminf_m F_12(e_(1/m)) >= h_tau-tau H.

Take the supremum over tau and use continuity/monotonicity of F_12 on
its branch. The stronger conclusion is

    liminf_m e_(1/m) >= e_*(H):=F_12^-1(D_*(H)) >0.              (3)

In particular m e_(1/m) tends to infinity. The source's finite liminf
bound rbar remains correct but is much weaker. Delta-dependent readouts
and a different stationary joining at each m do not evade(3): each joining
individually admits the decimation. A uniform finite upper bound on H also
suffices even for varying source flows; no construction of such sources
is claimed.

This is specifically a full stationary-history metric. It is NOT a claim
that finite-window or weak approximation is impossible. No bound uniform
over unbounded H, varying target p or independently retuned clocks follows.

## 3. Exact coarse control in the native g0 case

As a transparent native control, let g=0, p_i=1/12 and H<=1 nat per supplied
model-time unit. At tau=1/2, the off-diagonal entries are
(1-e^-1)/12<1/12. Their total probability is
(11/12)(1-e^-1)>11/24, since e>2. Also log12>2, since e<3.
Discarding the positive diagonal contribution gives

    h_(1/2)>11/12, h_(1/2)-H/2>5/12.

For d=1/100, h2(d)<=d log(1/d)+d and
F_12(d)<=d[log1100+1]<9/100<5/12. The log bound follows from
e>5/2 and (5/2)^8>1100. Therefore

    e_*(H)>1/100   for this g0/H<=1 control.                     (4)

This is intentionally coarse exact arithmetic, not a fitted threshold,
a bound for g=3 or measured device error. Positive-g numerical diagnostics
use their actual nonuniform p separately.

## 4. Compatibility controls

At one fixed delta an exact finite-entropy source with H=h_delta/delta
is not contradicted: h_(k delta)<=k h_delta, so all right sides in(1)
are nonpositive. Its entropy cost diverges upon refinement. Likewise a
clock-resolved exact CTMC path-space flow has infinite H and leaves the
finite-H class. A fresh random readout supplies additional entropy; its
cost cannot be attributed to a deterministic finite-H carrier for free.

Infinite-dimensionality, a singular observable or an exotic name alone
does not bypass a theorem that already permits arbitrary measurable
readouts on a standard probability space. The actual violated premise
must be identified. No source mechanism or physical law follows merely
from escaping this necessary obstruction.
