# RG-C03 — projected process and memory

For M=3, the sparse 1728-state CTMC is constructed directly from per-label heat-bath updates. Its coarse two-label observation has 144 values; the stationary conditional-averaged generator and the separately constructed marginal Gibbs sampler are both compared with actual projection. Rates for hidden states (0,0,0) and (0,0,1) differ, so strong lumpability **fails** at nonzero gain. The equilibrium distribution is shared by the two coarse generators, but finite-time transition distributions generally differ.

The scalar memory K(t)=⟨v,exp(t QLQ)v⟩_π, v=QLPf, f=1_{Y=(0,0)}, is evaluated for a reversible quadratic channel. K(0)=0.00995201523, derivative at 0 =-0.0180399818. A nonequilibrium hidden initial condition is recorded separately; it is not silently assumed stationary.
