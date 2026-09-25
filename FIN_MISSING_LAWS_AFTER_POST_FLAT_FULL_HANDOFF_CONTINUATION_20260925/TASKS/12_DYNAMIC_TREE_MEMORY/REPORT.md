# TREE-BULK-DYN-03 — dynamic Schur memory

Status: **PROOF_GRADE_CONDITIONAL**

## Main result
Giving internal tree nodes positive storage C_I turns static DtN reduction into a frequency-dependent boundary response.  Even one interior storage node generates a pole, so a frequency-independent coarse conductance cannot reproduce the dynamic response on an open frequency interval.

## Core formulas
\[
\Lambda(z)=L_{BB}-L_{BI}(L_{II}+zC_I)^{-1}L_{IB},
\]
\[
K_H(t)=L_{BI}e^{-C_I^{-1}L_{II}t}C_I^{-1}L_{IB}.
\]

## Evidence / reproduction
Algebraic Schur-complement calculation; this executes a concrete tree instance of the older REF-002 dynamic-subdivision program.

## Caveats
The storage law and clock remain added premises.

## Next question
Classify the pole/residue hierarchy and inverse reconstruction.
