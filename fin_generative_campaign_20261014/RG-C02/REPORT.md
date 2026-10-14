# RG-C02 — exact finite static reduction

The 4→3 elimination enumerates **all 1728 retained labelled states**, direct 4→3 sums of 12 and direct/sequential 4→2 sums of 144 removed states. FIN third contrast: **0.271270584589**. Pair-only approximation TV **0.0233588551**, residual oscillation **0.205388295**, universal bound tanh(osc/4) **0.0513019955**. The decomposition is uniform-reference ANOVA/Hoeffding; positive third term is an effective model term, not a source. The negative-gain check covers all three models.

Proof of TV bound: If r=F-Ftilde has range width ω, shift its constant so |r|≤ω/2. The likelihood ratio of normalized measures has bounded log-range ω; the maximal TV across positive reweightings with that log-range is tanh(ω/4). A constant residual has zero error. Applying any fixed Markov observation kernel contracts TV.

**Numerical status:** floating extrema are exhaustive over the finite alphabet, not machine-certified exact-kernel interval extrema; source-input enclosure is checked separately.
