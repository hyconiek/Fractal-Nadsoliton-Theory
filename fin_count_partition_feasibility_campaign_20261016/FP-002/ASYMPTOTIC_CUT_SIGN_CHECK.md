# Frozen phase cut: the macroscopic sign test is inconclusive, and why

Under the inherited MT critical limiting marginal and the frozen Gaussian radius cutoff, define `J_A(z)=int_{-infty}^z b(u)rho(u)p(u)du` and `J_B(z)=int_{-infty}^z b(u)rho(u)(1-p(u))du`, `b(u)=-u³/3`. For z<=0 every integrand is nonnegative, and strictly positive on a nonempty interval. For z>=0 use evenness of the phase density and oddness of b: `J_P(z)=-int_z^infty b(u)q_P(u)du>0`. These are strictly positive for every finite z.

This excludes only the **most elementary continuum parity-cut sign obstruction**. It does not prove finite-N `J_l^P>=0` for all levels, and still less the existence of a vertexwise nonnegative solution satisfying every microscopic incidence/cap constraint. In particular no Riemann-sum approximation near the moving finite-N phase boundary is claimed to yield a microscopic primal solution.

The phase-conditioned noise constants and the IF-feasible IBA lower bound are inherited from the user-provided BA-review seed (not available here as bytes): `P_A=1/2+E_rho exp(-Z²)/8`, `P_B=1-P_A`, `M=E_rho[Z²exp(-Z²)]>0`, `E[G|A]->2-M/(2P_A)`, `E[G|B]->2+M/(2P_B)`. Their numerical values in verify_fp.py are diagnostics only. This yields a potential nonzero bracket difference **only if** FP-002 is actually feasible. We do not activate FP-003 on this conditional statement.
