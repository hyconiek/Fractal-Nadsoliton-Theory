# RI-002 — The frozen 27-cell predictive experiment

**Scientific disposition:** finite-alphabet *risk convergence* and rigorous strict informativeness for the limiting experiment; numerical decimals are **uncertified finite-difference diagnostics**.

Fix `lambda=0,k=1,S=1`, previous integral observations at `t1=1/4,t2=1/2` and target `C=sign(Z_1)` with explicit third outcome `ZERO`. Define `V1=Y_{1/4}`, `V2=Y_{1/2}-Y_{1/4}`, and `Q(v)=LOW` for `v<=-1/2`, `MIDDLE` for `-1/2<v<=1/2`, `HIGH` otherwise. Record `R=(Q(V1),Q(V2))` has nine possible outcomes, all 27 target/record cells retained. The finite-N record uses the exact microscopic Y_N and target Z_N, with the same thresholds. No tuning or postselection.

**Limit of cells.** The admitted OBS joint J1 path convergence is to continuous `(Z,Y)` with `dZ=-Z^3/3 ds+sqrt(2)dB`, `dY=Z ds+sqrt(2)dW`, stationary quartic `Z0`, and independent B,W drivers. Each increment V_i has conditional law `N(A_i,2 Delta)` with `A_i=int_interval Z_s ds` and `Delta=1/4`; these conditional densities are continuous everywhere. Likewise Z1 has no atom at 0 (nondegenerate elliptic diffusion). Hence every quantizer-boundary and target-zero discontinuity is a null event. The continuous mapping theorem implies convergence of **all 27 cell probabilities**. On a fixed finite alphabet, coordinate convergence implies joint total variation convergence, even though the underlying path-TV and posterior convergence are not asserted. The total error retains the original iid-start path transfer and parent exit as o(1); no finite-N rate is claimed.

With `p_N(c,r)=P(C_N=c,R_N=r)`, the Bayes error for 0-1 loss and this **fixed nine-outcome record only** is

`r_N = 1 - sum_{r in {LOW,MIDDLE,HIGH}^2} max_{c in {MINUS,ZERO,PLUS}} p_N(c,r)`.

It converges to `r_* = 1 - sum_r max(p_*(-,r),p_*(+,r))` since `p_*(ZERO,r)=0`. Symmetry gives the no-record limit risk `1/2`.

**Exact limiting representation of the 27 cells.** Conditionally on the *signal path*, let `A1=int_0^.25 Z_s ds`, `A2=int_.25^.5 Z_s ds`, `sigma=sqrt(1/2)`, `H_q(a)=Phi((upper_q-a)/sigma)-Phi((lower_q-a)/sigma)` with standard infinite endpoints. Then

`p_*(c,(q1,q2))=E[ 1_{sign(Z1)=c} H_{q1}(A1) H_{q2}(A2)]`.

This is a complete 27-cell table *by exact integral definition*, not a value invented from N simulated samples. `p_*(ZERO,r)=0` and `p_*(-,(q1,q2))=p_*(+,(swap(q1),swap(q2)))` under global sign reversal. All cells are nonnegative and sum to one.

**Strict predictive advantage: rigorous qualitative theorem.** In the scalar SDE, solutions are strictly order-preserving with respect to the starting point for a common Brownian driver. Both `1_{Z1>0}` and `1_{V1>1/2}` are monotone functions of (a) the non-degenerate Z0 and (b) the underlying independent Gaussian driving increments (plus independent W observations). Scalar solution maps have nonnegative Malliavin derivative with respect to earlier signal noise (`sqrt(2) exp(int b'(Z)dt)>0`); the usual finite-time Gaussian association/Harris argument then gives nonnegative conditional covariance and strictly positive covariance from the varying initial Z0. Consequently `P(Z1>0,V1>1/2)>P(Z1>0)P(V1>1/2)`. The nine-cell table cannot be independent of the sign; its optimal risk is therefore **strictly less than 1/2**. This sign of advantage does not rely on uncertified numerical decimals. A quantitative nonzero advantage lower bound is not claimed.

**Deterministic numerical diagnostics** in `GRID_DIAGNOSTIC_COARSE.json` and `GRID_DIAGNOSTIC_FINE.json` solve a discretized reversible 1D diffusion in Z with Gaussian increments for Y, collapse increment bins at `1/4` and `1/2`, then propagate the sign to time1. These are NOT certified enclosures: truncation in Z/Y, birth-death spatial discretization, split-step time error, interpolation error and finite-grid atom at Z=0 are *not rigorously bounded*. The grid difference is a sensitivity diagnostic, not an error bar. Hence no decimal from these files is an accepted risk or 27-cell precision certificate.

Neither weak path convergence nor this fixed finite record proves mutual information convergence or sufficiency of all higher-resolution/parity histories.
