# RG-E03 fixed-ideal discrimination table

These numbers refer to the **ideal N2 difference law at g=3**, under trace-matched fixed matrices, not hardware data. The raw TV interval endpoints are independently computed under the mpmath.iv assumption.

| Fixed hypotheses | Ideal TV | Per-hypothesis max common TV error if neither has other error | Lower bound under *hypothetical* gain ±1% + instrument TV≤0.005 |
|---|---:|---:|---:|
| FIN / full Potts | 0.096998787035 | 0.048499394 | 0.074363614 |
| FIN / flat P7 | 0.011787191579 | 0.005893596 | 0.000000000 |
| Full Potts / flat P7 | 0.106340140736 | 0.053170070 | 0.083825732 |

**Readout nuisance rule:** if each allowed observed-law set is within TV radius r_i of ideal p_i, any pair has TV ≥ max(0,TV(p_i,p_j)-r_i-r_j). This follows from the triangle inequality and is a *guaranteed* minimax lower bound under that declared ball constraint, not an exact minimax optimum. Unknown arbitrary detector maps have lower bound zero. The ±1% and 0.005 figures are hypothetical design tolerances: NO measured calibration is present.

**Statistical error:** for 12 iid effective observations, a crude union-Hoeffding guarantee for TV empirical error≤eps at confidence≥1-alpha is n_eff≥(12²/(8 eps²)) log(24/alpha). At eps=0.002 and alpha=0.05 this gives 27,782,038 effective iid observations, only a conservative sufficient design number, not required/achieved sample size. Correlated trajectories demand separate category-specific mixing/ESS bounds or a predeclared valid block procedure; raw record counts are not iid samples.

An additional trace-preserving mode perturbation PERT_34_10P changes FIN's N2 law by TV≈0.014794589607. Arbitrarily smaller perturbations exist, so fixed-alternative separation never proves uniqueness in the unrestricted rank-seven class.
