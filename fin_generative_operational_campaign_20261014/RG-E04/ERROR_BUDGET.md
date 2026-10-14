# RG-E04 dynamic design error and kinetics scope

Exact model identities for all proposed q12 pair kinetics:
- A two-copy Gibbs equilibrium `pi(a,b)∝exp[(g/2)A_ab]` (constant diagonal).
- A single-copy heat-bath clock: at rate 1 resample a label from `softmax_j[(g/2)A_{j,other}]`.
- Competing proposals to each of the 11 different labels, with Metropolis acceptance `min(1,exp(delta))` and attempt rate `11/12` per copy; Barker acceptance `1/(1+exp(-delta))` and attempt rate `22/12` per copy.
- At g=0 the two adjusted competing generators are exactly the same as heat-bath: per offdiagonal single-site transition rate `1/12`; the nonconstant single-label Fourier correlation is `exp(-t)` in dimensionless time. This is a **defined model calibration**, not evidence that any device retains the same physical clock when interaction is switched on.

The prospective **N2** prediction fixes uniform preparation at g=0, switch to g=3 at model t=0, a noninvasive difference readout at t1=0.4, and another at t2=1.0, with outcomes (0..11, INVALID). The full 13x13 joint is computed from the 144-state parent; no M3 ROM is used. All nominal invalid mass is zero because the *ideal* instrument is defined that way.

**Full physical uncertainty budget is OPEN**, not estimated from simulations:
1. actual finite channel compilation/gain and interaction drift;
2. detector confusion, latency, back-action and invalid/no-hot/multi-hot mass;
3. g=0 calibration time conversion and evidence of rate constancy under g=3;
4. preparation fidelity, nonstationarity, record dependence/ESS and adaptation effects;
5. numerical expm_multiply floating roundoff (independently checked with dense expm; no outward matrix-function certificate);
6. alternative reversible kinetics beyond the three *declared* laws.

No physical PASS or uniform error interval is assigned. The previous M3 heat-bath hidden spectral bound I≤H≤3I is **not** transferred to N2 or alternative kinetic laws. Full multinomial/time-series likelihood analysis requires data with timestamps, invalid records and independently frozen nuisance models.
