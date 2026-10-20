# AH-002 — fully quantitative finite-prefix recovery from equality names

Date: 2026-10-20. Activated after AH-001's positive theorem.

## Assumptions and distance

Use the identical 12-state jump-chain family and first-appearance encoding specified in AH-001. Assume a **declared lower bound** `J_ij>=a>0` for every `i!=j`, `a<=1/11`. The lower bound is analysis-class information, not a sourced microscopic FIN law. Let

`d_iso(A,B)=min_{P permutation} ||A-P B P^T||_max`.

For two supplied hypotheses `J,K`, assume `d_iso(J,K)>=delta>0` (`delta<=1`), and choose any target `eta in (0,1)`. Put `b=10 a^2` and

```
m = ceil( max{ 2/b,
               8/b^2 * log(24/eta),
               384/(delta^2*b^2) * log(528/eta) } )
L = 2*m+2   observed states / (2*m+1) observed transitions.
```

**Theorem AH2.** There exists a deterministic decision rule using only the first-appearance equality prefix of length `L` whose probability of choosing the wrong hypothesis is at most `eta` under **each** of the two models. More generally, the same construction yields a high-probability `d_iso` reconstruction bound for an unknown admitted `J`.

### Proof: coverage and row denominators

Write `n=L-1=2m+1`; count departures from all states at times `t=0,...,2m`. For every original vertex `i` and each `k=1,...,m`, the two-step lower bound yields

`E[1_{Y_{2k}=i} | Y_0,...,Y_{2k-2}] = J^2(Y_{2k-2},i) >= b`.

After subtracting conditional expectations, these are martingale differences of absolute value at most 1 for the even-time filtration. Azuma-Hoeffding gives

`Pr(sum_{k=1}^m 1_{Y_{2k}=i} < b*m/2) <= exp(-b^2*m/8)`.

For all twelve original vertices simultaneously, except on a set of probability at most `12 exp(-b^2 m/8)`, the departure counts `N_i>=b*m/2`. The choice `m>=2/b` makes all 12 counts positive and guarantees coverage in the observed prefix. This step does **not** treat visits as independent Bernoulli variables.

### Proof: edge-count concentration

For each ordered `i!=j`, let

`M_ij=sum_{t=0}^{n-1} (1_{Y_t=i,Y_{t+1}=j}-J_ij 1_{Y_t=i})`.

These are sums of martingale differences of magnitude <=1 (conditional on the previously seen trajectory). Azuma-Hoeffding yields, for any `u>0`,

`Pr(|M_ij|>=u) <= 2 exp(-u^2/(2n))`.

Choose `u=delta*b*m/8`. Since `n=2m+1<=3m`,

`Pr(|M_ij|>=delta*b*m/8) <= 2 exp(-delta^2*b^2*m/384)`.

There are 132 ordered off-diagonal pairs; taking their union costs at most

`264 exp(-delta^2*b^2*m/384)`.

When all occupation and edge events hold, with `N_i>=b*m/2`, the empirical transition matrix `\hat J_ij=N_ij/N_i` satisfies

`|\hat J_ij-J_ij| <= (delta*b*m/8)/(b*m/2)=delta/4`.

Both empirical and model diagonals equal zero, so this bound holds for every entry. The first-appearance encoding is simply an unknown random simultaneous permutation of every visited row and column, not an additional source of statistical error. Thus

`d_iso(\hat J, J_true)<=delta/4`.

By the chosen lower bound on `m`, each union term is <=eta/2. Hence the claimed high-probability reconstruction event has probability at least `1-eta`.

### Proof: deterministic classifier

Construct `\hat J` solely from anonymous transition/occupation counts; if some row was never observed, output either fixed hypothesis arbitrarily. Otherwise return the hypothesis minimizing `d_iso(\hat J,H)` over `H in {J,K}`. On the high-probability event, the true model has distance <=delta/4. The other has distance >=delta-delta/4=3delta/4 by the triangle inequality for the quotient metric. Therefore the rule chooses correctly, under either model, with error probability at most `eta`. QED.

### Limitations

- This bound is **sufficient, not optimal**. A demonstration with `a=1/20`, `delta=1/10`, `eta=1/20` exceeds a billion sampled *states*. This is not a feasible sample-budget recommendation and is not meant to model physical measurements.
- The best permutation in `d_iso` defines the decision rule mathematically; the argument does **not** enumerate 12! cases. Finite group size guarantees existence of the minimum. No practical algorithmic complexity guarantee is claimed.
- No independence between successive jumps was assumed. All randomness and correlation are those of the one supplied homogeneous Markov chain.
- The quantitative bound applies only when one has a valid **known** `a` and fixed quantitative separation `delta`. It does not infer these from FIN's information ontology.
