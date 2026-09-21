# R7N-044 K20 pivot decision

The K16 checkpoint path was externally restored to the earlier 114000-cell snapshot after a later 200000-cell run. The textual logs for 119000--200000 remain, but the complete later partition tree is unavailable and therefore is **not** accepted as a proof object.

The campaign resumes from the intact 114000-cell K16 proof checkpoint. This is not a repetition of the same cover architecture: K20 reduces the rigorous uniform full-gradient remainder bound from approximately 2.38e-8 to approximately 1.99e-10 while retaining only 45 resonances. The next bounded pass first reclassifies the 5382 intact K16 residual cells without subdivision, then permits root-aware K20 subdivision only for the remaining cells.

This pivot is justified by the >100x tighter certified remainder and avoids reconstructing a lost later K16 tree merely to recover diagnostic progress.
