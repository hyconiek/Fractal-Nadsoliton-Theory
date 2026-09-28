# FIN 332 — PREPARATION-CONTROLLER-WITHOUT-BASIN-ORACLE
## Unrestricted biased heat-bath can prepare the localized phase without a hard J=0 wall

Date: 2026-09-28

Status: **PASS as an operational soft-confinement mechanism for the already-opened finite-N range, with explicit escape outcomes; low-N escape is non-negligible and must not be hidden.**

## 1. Controller

Use the unrestricted leave-one-out heat-bath with field on label 0:

`q_j^prep(n,i) = softmax_j[(g/N) A(n-e_i) + theta e_0]`.

No knowledge of the basin label `J=0` is used by the controller.

The invariant distribution is exactly the unrestricted biased Gibbs measure

`pi_theta(n) proportional to pi_N(n) exp(theta n0)`.

The hard-wall target used previously is simply this distribution conditioned on `J=0`.

Therefore the exact TV distance between unrestricted equilibrium and the hard-wall target is

`1 - pi_theta(J=0)`.

## 2. Static soft confinement at theta=2

At `theta=2`, the equilibrium mass in the desired phase is:

- N=3: 95.3170%;
- N=4: 99.1376%;
- N=5: 99.8868%;
- N=6: 99.9808%;
- N=7: 99.99756%;
- N=8: 99.99957%;
- N=9: 99.999946%;
- N=10: 99.999990%.

Thus for N>=5 the field itself is already an effective soft confinement. No postselection is required if other phases and `unlocalized` are retained as outcomes.

## 3. Leakage from the phase

The deep seed has exactly zero instantaneous one-step exit rate from J=0 for every tested N=3..10.

The mean exit flux evaluated under the conditional target falls rapidly:

- N=3: `1.965e-2`;
- N=4: `4.976e-3`;
- N=5: `7.151e-4`;
- N=6: `1.793e-4`;
- N=7: `2.342e-5`;
- N=8: `6.013e-6`;
- N=9: `7.646e-7`;
- N=10: `2.000e-7`.

The large maximum boundary rates are irrelevant for seed preparation because the target assigns extremely small mass to those boundary states.

## 4. Exact no-exit calculation from the deep seed

A killed subgenerator was built on J=0. Its diagonal includes all exit rates, so

`S(t)=sum p_killed(t)`

is the exact probability that the unrestricted controller has **never left the basin** before time t.

At `Tprep=4`:

| N | escape by t=4 | TV(survivor,target) | conservative total bound |
|---:|---:|---:|---:|
| 3 | 4.4910% | 1.4166% | 5.8439% |
| 4 | 1.2418% | 0.5146% | 1.7500% |
| 5 | 0.1436% | 0.1754% | 0.3188% |
| 6 | 0.03760% | 0.1119% | 0.1494% |
| 7 | 0.00417% | 0.09072% | 0.09489% |
| 8 | 0.00110% | 0.08727% | 0.08837% |
| 9 | 0.000123% | 0.08715% | 0.08727% |
| 10 | 0.0000329% | 0.08867% | 0.08871% |

The last column is

`P(exit by t) + P(survive to t) * TV(survivor,target)`.

It is a valid conservative bound if every escaped trajectory is treated adversarially.

## 5. Meaning

For the finite-N regime where the effective process is already clean (roughly N>=5, especially N>=7), the hard basin oracle is not needed operationally.

A simple physical-style controller exists:

`deep seed -> turn on field theta -> relax -> turn field off -> retain every outcome`.

The field does not need to know which microstates belong to J=0.

At N=3,4 the escape budget is still material. Those cases must be reported explicitly rather than postselected away.

## 6. Verdict

**PASS for controller realizability without a hard basin wall in the studied larger-N regime.**

This does not source a laboratory actuator or physical meaning for `theta`; it only supplies a complete operational stochastic protocol inside the declared finite model.
