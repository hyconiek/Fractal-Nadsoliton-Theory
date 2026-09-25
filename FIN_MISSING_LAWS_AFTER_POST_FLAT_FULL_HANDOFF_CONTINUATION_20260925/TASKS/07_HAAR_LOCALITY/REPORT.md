# LOCALITY-BRIDGE-01 — hierarchy gradient and projective consistency

Status: **PROOF_GRADE_CONDITIONAL**

## Main result
Given a hierarchy H, normalized Haar contrasts define a hierarchical gradient `B_H`.  With positive scale weights `W`, the boundary operator is `L_H=B_H^* W B_H`.  Coarse fields copied identically to children generate no new high-frequency Haar modes, giving exact projective consistency under refinement.

## Core formulas
\[
E_H(f)=\|W^{1/2}B_Hf\|^2,\qquad L_H=B_H^*WB_H,
\]
\[
I^*L_{L+1}I=L_L.
\]

## Evidence / reproduction
Algebraic; regular-tree numerical checks were at machine precision.

## Caveats
This is locality in hierarchy-mode space, not yet ordinary nearest-neighbor physical space.

## Next question
Embed the boundary form into an explicitly local bulk network.
