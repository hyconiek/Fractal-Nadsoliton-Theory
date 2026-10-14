# TH-004 auxiliary-vs-actual dictionary

| Object | Meaning | Path status |
|---|---|---|
| `theta_ret=g X^T p_N` | observable of the leave-one-out jump chain | actual microscopic collective observable |
| `theta_aux` | Hubbard–Stratonovich latent Gaussian | equilibrium representation only |
| `-g grad Phi` | leading drift of theta_ret | derived from supplied heat bath |
| `g^2 X^T D(p)X/N` | leading covariance rate of theta_ret | derived from jumps |
| `2L/N` in MP7-029 | covariance of a separately supplied Langevin model | not transferable to actual labels |
| `sqrt(2g_f)` soft amplitude | fold projection of actual jump covariance | derived identity, source-certificate conditioned |
