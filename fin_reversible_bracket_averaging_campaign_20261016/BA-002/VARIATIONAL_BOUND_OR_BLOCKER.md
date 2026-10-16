# BA-002 — Variational / Poisson estimate and precisely unpaid q12 bound

## Exact finite reversible statement (rederived)

For a finite possibly reducible stationary reversible `L`, let `r=G-E_pi G` and `h=P_kerL r`. The componentwise mean-zero Poisson equation `-L u=r-h` has a unique solution after choosing zero mean **separately on each communicating component**. Its cost is

`R=<r-h,u>_pi = <u,-L u>_pi = sup_f [2<r-h,f>_pi - D(f,f)]`, with `D(f,f) = sum_{x<y} c_xy(f_y-f_x)^2`.

All quantities are finite for each fixed N even when the generator is reducible; there is no inverse on `h`. For stationary time `[0,T]`, the forward Dynkin martingale `M^+_t=u(X_t)-u(X_0)+integral_0^t(r-h)(X_s)ds` and the *reversed*, backwards martingale `M^-_s` satisfy for fixed `T`

`2 integral_0^t(r-h)(X_s)ds = M^+_t + M^-_T - M^-_(T-t)`.

Because `E sup|M^+|²` and `E sup|M^-|²` are each at most `8 T R` by Doob and `E Gamma(u,u)=2R`,

`E sup_(t<=T)|integral_0^t(r-h)(X_s)ds|² <= 20 T R`.

The invariant part is constant along stationary paths: `integral_0^t h(X_s)ds = t h(X_0)`. Split the parts with `|a+b|²<=2a²+2b²`, set `T=sqrt(N)S` and divide by `sqrt(N)` for scaled integrals to obtain

`E_pi sup_(s<=S)|integral_0^s r(X_(sqrt(N)v))dv|² <= 2 S² ||h||²_pi + 40 S R/sqrt(N)`.

Thus `||h||_pi->0` and `R=o(sqrt(N))` suffice for IBA, exactly as the user-pasted admitted review formula. The displayed constant is nonoptimal and not a convergence rate from the frozen class alone.

## New sector decomposition

Let `P_m` be equilibrium conditional expectation onto parity levels, `r_parallel=P_m r` and `r_perp=r-r_parallel`. BA-001 gives the exact finite-N formula `r_parallel=P_m(G_native-E_pi G_native)` and hence `||r_parallel||_L1(pi)->0`. Because `|r|<=8C` uniformly, `||r_parallel||_2->0` also. This settles **all static parity-visible noise errors** without NC.

It does **not** pay `R`, because the full Dirichlet-form dual `D(f,f)` ranges over functions of the entire count state, and `r_perp` can correlate with functions that distinguish slowly communicating configurations within the same parity level. The identity `D_tilde(f(m),f(m))=D_native(f(m),f(m))` provides no inequality of the form `D_tilde(f,f)>=const D_native(f,f)` for arbitrary `f`. In particular one cannot apply the known native HB mixing bound to a diagnostic chain with arbitrarily small activity.

## Rigorous rational caution: near-disconnected six-node stress

The verifier uses a fixed stationary six-vertex cycle with conductances

`c(tau)=(tau,tau,1/5-tau,1/5-tau,tau,1/5-tau), 0<tau<=1/10`.

At every tau it has the SAME invariant measure and the SAME drift of the displayed three-level parity variable. At tau=0 the positive-edge components are `{0,5}`, `{1}`, `{2,3,4}` and the invariant-component projection of the centered bracket is nonzero. For tau>0 it is irreducible (so `h=0`), but the exact rational Poisson cost satisfies

`lim_(tau downarrow 0) tau R(tau) = 1064/126075 > 0`.

This limit follows from finite-dimensional singular perturbation: the only `1/tau` part of `(-L)^(-1)` is on the three-component equilibrium-mean subspace. The reduced triangle has edge conductances `tau` and masses `(15,3,23)/41`; its Poisson equation `K u=diag(masses) h` with `K=[[2,-1,-1],[-1,2,-1],[-1,-1,2]]` has energy `h^T diag(masses) u =1064/126075`. All within-component inverses stay bounded because their positive edges tend to `1/5`. This computes the exact leading coefficient without relying on fitting to numerical points.

The fixture **does not** satisfy the q12/N critical sequence and cannot itself disprove universal q12 IBA. It does prove that drift, equilibrium, one-step support, bounded rates and parity Dirichlet information alone do not yield an abstract uniform Poisson bound on unrestricted finite graphs.

## Status of the native variational obligation

No estimate uniform over the **actual** all-N q12/A7 admissible conductance polytope has been obtained for either `||h_N||->0` or `R_N=o(sqrt(N))`. A new class-specific geometric/flux coercivity theorem, or a valid alternate maximal additive-functional argument, remains required. This is an explicit `OPEN_VARIATIONAL_Q12_CONDUCTANCE_BOUND`, not a conclusion that universal IBA is false.
