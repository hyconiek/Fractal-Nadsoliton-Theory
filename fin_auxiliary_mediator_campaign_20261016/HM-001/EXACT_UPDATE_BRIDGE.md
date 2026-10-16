# HM-001 exact update bridge

The normalizing Gaussian integral for candidate j yields `exp[g/(2N)||X^T(r+e_j)||²]`. The j-dependent part is `(g/N)(Ar)_j+(g/(2N))A_jj`; constant diagonal cancels. Consequently an ephemeral ν_r draw followed by categorical conditional output gives exactly the original LOO refresh law, including null refreshes. Any attempt to integrate using ρ_n instead changes the kernel. See HM_THEOREM_PROOFS.md §HM-001 for the full derivation.
