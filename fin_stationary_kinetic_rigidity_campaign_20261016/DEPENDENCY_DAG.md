# Dependency and evidence DAG

```
q12/A7 supplied exact Gibbs and reflected native generator [predecessor]
    + UF-002 exact scalar drift remainder [predecessor]
    + MT-003 stationary quartic/Gaussian moments [predecessor]
    + MT boundary geometry/exponential mass [predecessor]
    + diagnostic reversible count CTMC, CN rate cap, exact ALL-STATE drift match [KR class]
       -> stationary identity E Gtilde = E Gnative -> 2 [KR-002]
       -> vanishing stationary scaled drift error [KR-002]
    + NC [INDEPENDENT SUFFICIENT ASSUMPTION]
       -> integrated bracket 2t, joint martingale convergence [KR-002]
       -> stationary critical parity diffusion [CONDITIONAL THEOREM]

without NC:
    exact stationary generator balance + UF/MT -> E[h(Z_N)(G_N-2)] ->0 [KR-003]
    but full trajectory bracket -> 2t [OPEN IBA]
```

Intrinsic count statistics only; no kinetic-source bridge activation.
