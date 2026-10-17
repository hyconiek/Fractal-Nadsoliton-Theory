# UE-004: quantitative nonautonomous coarse-map obstruction

## Frozen coarse interface

Take the declared local Cauchy phase space `X=(Phi,Phi_t,Psi,Psi_t)` with the needed Sobolev regularity on the 1D line. Define the **fixed** coarse map

`C(X)=(Phi,Phi_t)` in `H^s × H^{s-1}`.

The hidden coordinates `Psi,Psi_t` are discarded. The reduced predictor would be any deterministic autonomous evolution of `C(X)` alone. This is not a fitted effective model. For the quantitative obstruction observe the `Phi` component in `L²(R)`; require sufficiently regular bounded solutions up to a common `T` and a uniform third-time-derivative remainder bound `||R_i(t)||_L2 <= C t^3` for `0<=t<=T`.

## Equal coarse state / different acceleration theorem

Choose two valid initial full states `X_1,X_2` with **identical** `(Phi_0,Phi_t(0))` and hidden data yielding

`a=|| Phi_tt^1(0)-Phi_tt^2(0) ||_L2 >0`.

The evolution yields, by Taylor and remainder bound,

`||Phi^1(t)-Phi^2(t)||_L2 >= (a/2)t² - 2C t³`.

For `0<t<= min(T,a/(8C))` (with the natural interpretation if `C=0`), this gives

`||Phi^1(t)-Phi^2(t)||_L2 >= (a/4)t²`.

Any deterministic autonomous coarse predictor `P_t(C(X))`, seeing identical coarse initial data, must therefore have worst-case error at least

`max_i ||Phi^i(t)-P_t(C(X_i))_Phi||_L2 >= (a/8)t²`.

This is a rigorous **finite-time, quantitative obstruction** to exact autonomous closure on a state class containing both preparations, under the stated smooth-solution bounds. It does not exclude a reduction with memory, retained hidden data/momenta, a restricted preparation manifold or a controlled error exceeding this lower bound.

## Exact toy control showing the obstruction is nonvacuous

In the **illustrative** `V_int=Phi²Psi²` scalar action, fix `Phi_0=tanh x`, `Phi_t(0)=0`, and choose the hidden profiles `Psi_0^1=0`, `Psi_0^2=epsilon sech x`, both with `Psi_t(0)=0`. Because `Phi_tt=Phi_xx-U'(Phi)-2Phi Psi²`, acceleration difference is exactly

`Delta a(x)=-2 epsilon² tanh x sech² x`, with `||Delta a||_L2=4 epsilon²/sqrt(15)>0`.

The same initial `Phi` and `Phi_t` yield different `Phi` futures. The exact integral is checked in `REPLAY/verify.py`. This is not a FIN-specific coefficient test. To upgrade it, bind the actual action and exhibit a nonzero `partial_Phi V_int` difference in its admissible full complex sector.

## Architecture / scaling stop

- Retaining `(Phi,Phi_t,Psi,Psi_t)` restores trivial exact autonomy; that is no coarse reduction.
- A projection onto `Phi` only has no universal Markov/second-order closed law when the hidden force depends on `Psi`.
- A history/memory closure must define the retained time window and prove an approximation estimate, not infer it from local Taylor algebra.
- Fixed kink width, exponent `eta=1.8`, twelve labels and q12 finite-size powers do **not** provide a family of nested spatial scales. A spatial fractal claim requires explicit refinement maps `R_{n+1->n}`, energy/observable consistency, admissible topology and quantified errors across scales. None is supplied here.

**Decision:** `QUANTIFIED_COARSE_AUTONOMY_OBSTRUCTION_IN_CONDITIONAL_CLASS; NO_REFINEMENT_SOURCE`.
