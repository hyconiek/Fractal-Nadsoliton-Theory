# Results and interpretation

The refined sampled maximum energy gain is 7.650175164 at t=116.5. It falls to 7.113042992 at t=200. Thus the previously endpoint-limited growth is followed by a decline within the extended observation window. This supports a transient-amplification interpretation over this window, not a theorem of asymptotic saturation or global stability.

| t | Optimized gain | Fixed position seed | Fixed velocity seed |
|---|---|---|---|
| 0 | 1.000000000 | 1.000000000 | 1.000000000 |
| 50 | 5.066685188 | 3.990313806 | 5.040409694 |
| 100 | 7.522447255 | 5.575869596 | 7.522187129 |
| 150 | 7.387348196 | 5.548858880 | 7.387251847 |
| 200 | 7.113042992 | 5.331611443 | 7.112851622 |

Both fixed seeds grow and later decline, so the observed pattern is not solely caused by changing the optimizing initial phase.

Maximum absolute gain-trace difference under dt halving: 5.007926264788409e-06. Historical [0,100] overlap difference at equal dx/dt: 1.4210854715202004e-14. Static energy drift: 8.371858761790918e-12. These are numerical checks, not rigorous error bounds. No spatial refinement of the extended horizon was executed.

Next bounded step: compute the signed reference-energy work budget on this same background. For each fixed perturbation, dE_ref/dt=integral[(V_ref-V(t)) z z_t]. Separate the phi-background contribution from 3 psi^2+5 psi^4, integrate at RK4 stages, and compare to the independently measured energy change. This tests the mechanism of growth and decline without launching another unbounded horizon scan. No exponential-instability or global-stability claim is warranted.
