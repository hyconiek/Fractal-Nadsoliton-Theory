# FIN 330 — PREPARATION-GENERATOR-REALIZABILITY
## From a target preparation distribution to a finite-time preparation process

Date: 2026-09-28

Status: **PASS for exact unconditioned biased Gibbs stationarity; PARTIAL for basin-conditioned realization; PASS for seed-conditioned operational mixing over N=3..10 at theta=2.**

## 1. Biased leave-one-out generator

Let

`theta = kappa/N`.

For count state `n`, remove one copy from label `i`, obtaining `m=n-e_i`, and reinsert it at label `j` with

`q_j^prep(m) = softmax_j[(g/N)(A m)_j + theta 1_{j=0}]`.

The count-state transition rate is

`Q(n,n-e_i+e_j) = n_i q_j^prep(n-e_i)`.

## 2. Exact detailed balance for the full biased Gibbs law

The candidate stationary distribution is

`pi_{N,kappa}(n)`

` proportional to N!/prod_i n_i!`

` * exp[(g/(2N)) n^T A n + (kappa/N)n0]`.

For neighboring states `n=m+e_i` and `n'=m+e_j`, the ratio of multinomial factors gives `n_i/(n_j+1)`. Because the FIN matrix is circulant, all diagonal entries `A_ii` are equal, so the quadratic-energy difference is exactly

`(g/N)[(A m)_j-(A m)_i]`.

The pinning difference is `theta(1_{j=0}-1_{i=0})`, exactly matching the softmax ratio. Hence

`pi(n) Q(n,n') = pi(n') Q(n',n)`.

Therefore the **unrestricted** biased leave-one-out heat-bath chain realizes the full biased Gibbs distribution exactly.

## 3. Why basin conditioning is not automatic

The target used in the recent FIN preparation contract is not the full biased Gibbs law. It is the conditional law

`mu_{N,kappa} = pi_{N,kappa}(. | J=0)`.

A simple way to preserve it is to reflect every proposed move that exits the declared `J=0` basin. Detailed balance then remains true on every retained edge.

However, the resulting basin graph is not always connected. Number of connected components:

- N=3: 1;
- N=4: 2;
- N=5: 2;
- N=6: 9;
- N=7: 16;
- N=8: 34;
- N=9: 41;
- N=10: 112.

Thus a reflecting chain started from one state does **not** in general converge to the full basin-conditioned distribution. It converges to that distribution further conditioned on its connected component.

This is a genuine no-go for the naive statement “add the field and reflect at the basin boundary, therefore the exact conditional preparation is realized.”

## 4. The deep-seed component nevertheless carries almost all target mass

Let the operational initial state be the deep seed `n=(N,0,...,0)` and let `C_seed` be its connected component inside the reflecting basin graph.

At `kappa=0`, the mass outside `C_seed` is already small for N>=5; for example:

- N=6: `5.1361e-4`;
- N=7: `9.8256e-5`;
- N=8: `7.7738e-5`;
- N=9: `1.7176e-5`;
- N=10: `1.2888e-5`.

For the cross-N-stable control `theta=2` (`kappa=2N`), it becomes tiny:

- N=4: `1.4665e-6`;
- N=5: `1.7634e-9`;
- N=6: `4.4933e-9`;
- N=7: `5.9060e-10`;
- N=8: `1.2649e-11`;
- N=9: `1.3430e-12`;
- N=10: `4.8628e-14`.

The TV distance between the full basin-conditioned target and the target additionally conditioned on `C_seed` is exactly this omitted component mass.

Thus reflecting dynamics + deep-seed initialization is an extremely accurate realization of the declared conditional target in the tested positive-bias regime, even though it is not mathematically identical to it.

## 5. Spectral gaps are not the operational mixing time

At `theta=2`, the reversible reflecting generator on `C_seed` has the following spectral gaps:

- N=3: 0.95196;
- N=4: 0.91672;
- N=5: 0.35601;
- N=6: 0.92472;
- N=7: 0.37018;
- N=8: 0.66736;
- N=9: **0.0372457**;
- N=10: 0.40242.

The N=9 value gives a very slow worst-state relaxation time `1/gap ~=26.85`. Its two leading nontrivial eigenvalues are nearly degenerate.

A direct eigenvector check shows, however, that the deep seed has essentially zero overlap with this anomalous slow pair: approximately

`-3.5e-8` and `6e-13`

for the two normalized slow eigenfunctions.

Therefore the worst-state gap is not the relevant preparation time for the declared operational initial condition.

## 6. Direct seed-to-target mixing

The exact finite-state semigroup was propagated from the deep seed for every N=3..10 at `theta=2`.

Grid-resolved results:

- `TV < 1%` by `t ~=1.5` for N=6..10 and by `t ~=1.75` for N=3..5;
- `TV < 0.1%` by `t ~=4.0` for N=6..10 and by `t ~=4.25` for N=3..5.

At N=9 specifically:

- t=1: TV `0.01571`;
- t=2: TV `0.00596`;
- t=4: TV `0.000871`;
- t=8: TV `1.88e-5`.

The N=9 worst-state spectral anomaly therefore does **not** prevent rapid preparation from the deep seed.

This gives a useful operational distinction:

- **uniform worst-state mixing theorem:** not established and clearly poor at N=9;
- **deep-seed preparation protocol:** numerically fast and stable over N=3..10 for `theta=2`.

## 7. Remaining physical/operational gap

The reflecting boundary uses the already declared basin label `J=0`. Implementing an exact reflecting wall may itself require an oracle-like controller that recognizes every attempted basin exit. This is not yet sourced by FIN.

More physical alternatives remain to be tested:

1. unrestricted biased preparation followed by an explicit localization readout/postselection;
2. a soft confining controller rather than a hard basin wall;
3. a feedback rule based only on available observables, with an explicit information/control cost.

If postselection is used, an `unlocalized/rejected` outcome should be kept explicitly in the observed process or its probability paid in the error budget. Conditioning it away is not a free TV contraction.

## 8. Main conclusion

The preparation problem now has a concrete chain:

`control theta -> exact biased heat-bath -> (approximate) basin restriction -> finite mixing time -> prepared microscopic law -> observed process`.

This is materially closer to an operational statistical-mechanics theory than an abstract fitted preparation vector. It still does not source a physical actuator, clock unit, or measurement apparatus.
