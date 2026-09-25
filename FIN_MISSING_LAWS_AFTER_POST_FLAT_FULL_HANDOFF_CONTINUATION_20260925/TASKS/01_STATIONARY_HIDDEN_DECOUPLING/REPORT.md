# EDGEWORTH-MEM-03 — stationary hidden OU decoupling

Status: **PROOF_GRADE_CONDITIONAL**

## Main result
At any stationary declared heat-bath state `p*=q*`, write visible/hidden Gaussian fluctuations `(x,y)` and define `K=D_yx D_xx^{-1}` and `z=y-Kx`.  The algebraic identity `C=K(A+I)` makes `z` an independent OU residual with drift `-z`.  Therefore leading Gaussian visible dynamics is Markov and the first genuine equilibrium hidden-memory feedback appears only after two hidden couplings, i.e. at order `1/N`.

## Core formulas
\[
F=X^TSX,\ H=Y^TSX,\ G=Y^TSY,\quad S=\operatorname{diag}p-pp^T,
\]
\[
A=gF-I,\ C=gH,\ D_{xx}=2F,\ D_{yx}=2H,\quad K=HF^{-1},
\]
\[
z=y-Kx,\qquad dz=-z\,dt+dW_z,\qquad \operatorname{Cov}(dW_z,dW_x)=0.
\]

## Evidence / reproduction
Recomputed directly from strict X7 and orthonormal k=1,2 Fourier hidden modes in `REPLAYS/replay_strict_memory.py`.

## Caveats
Heat-bath is declared; outside stationarity the exact decoupling need not hold and `N^-1/2` hidden dependence can reappear.

## Next question
Derive the complete reduced `1/N` Edgeworth operator, including local and memory terms.
