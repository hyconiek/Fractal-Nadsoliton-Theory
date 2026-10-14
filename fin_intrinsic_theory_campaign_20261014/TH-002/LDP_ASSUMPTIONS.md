# TH-002 LDP assumptions and statuses

| Kraaij Theorem 2.3 premise | FIN status |
|---|---|
| finite q-state mean-field empirical chain | PASS, q=12 |
| continuous limiting `v_ij(p)` | PASS, analytic softmax |
| positive when source occupation positive | PASS |
| boundary monotone factorization | PASS: `v=p_i*s_j` |
| uniform finite-N rate convergence | PASS, O(1/N) for fixed/compact g |
| initial-law good LDP | CONDITIONAL generally; deterministic and equilibrium Gibbs classes paid |
| topology | Skorokhod path space on compact simplex, as in theorem |

The theorem proves a path LDP for the declared mathematical chain. It says nothing about physical time, a natural source of A7, or experimental data.
