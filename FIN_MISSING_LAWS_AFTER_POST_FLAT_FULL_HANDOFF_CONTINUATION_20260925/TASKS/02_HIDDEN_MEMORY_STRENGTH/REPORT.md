# MEMORY-KERNEL-04/05 — residual covariance and memory strength

Status: **RECOMPUTED_NUMERIC_PLUS_FORMULA**

## Main result
The residual hidden covariance is the Schur complement `Sigma_z=G-HF^{-1}H^T`.  The hidden modulation of visible diffusion is `delta D=sum_a z_a T_a`, `T_a=X^T diag(Y_a) X`.  A basis-independent strength `M=E||delta D||_F^2` decreases strongly from uniform to localized phase.

## Core formulas
\[
\Sigma_z=G-HF^{-1}H^T,\qquad
\mathcal M=\sum_{ab}(\Sigma_z)_{ab}\langle T_a,T_bangle_F.
\]
Stationary autocorrelation: \[
\langle\delta D(t),\delta D(0)angle_F=\mathcal M e^{-|t|}/N+\cdots.
\]

## Evidence / reproduction
Recomputed in `REPLAYS/replay_strict_memory.py` from strict model inputs.

## Caveats
The unit decay rate is in the declared heat-bath clock normalization.

## Next question
Connect this channel to full quartic/Edgeworth observables.
