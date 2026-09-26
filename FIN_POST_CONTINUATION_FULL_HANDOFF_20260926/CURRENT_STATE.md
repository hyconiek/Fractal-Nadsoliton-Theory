# Current state after all post-checkpoint research

## A. Finite-N heat-bath / hidden-memory lane

**Closed within the declared model:**
- KURTOSIS-THEOREM-11: pure k=5 coefficient `lambda5^2(g lambda5-6)/144` with exact Z12 derivation; k=3,4,6 pure self-quartics are null controls.
- GENERAL-QUARTIC-THEOREM-13: for any retained real probe `phi`, with `h=P_H(phi^2)` and `k=P_H(phi A7 phi)`, the leading defect is
  `C(phi)=-12<h,h>_u+2g<h,k>_u`; the `g^2` tensor vanishes.
- HIDDEN-MULTIPLICATION-RANK-14: `Sym^2(V7)->H4` has exact rank 4; all four discarded k=1,2 real modes are quadratically reachable.
- EDGEWORTH-STATIONARY-15: stationary reduced generator is closed through O(1/N) in the declared lane; the leading Gaussian hidden residual decouples and the first genuine reduced equilibrium feedback occurs at O(1/N).
- FINITE-N-REFRESH-INVARIANT-16: exact invariant/reversibility statement for the declared refresh chain.
- FIBER-EDGEWORTH-BIAS-17: finite-N conditional fiber bias relative to ME7 is quantified.
- FINITE-N-TARGET-CORRECTION-THEOREM-21 and GENERAL-RATE-CORRECTION-23: leading quartic convention sensitivity depends only on a four-dimensional hidden net first-jump drift; zero net drift and retained-only drift are invisible.

**Probe design:**
- Mixed retained probes can be much stronger than pure k=5.
- QUARTIC-PROBE-INTERVAL-CERTIFICATE-23 rigorously bounds the global maximum value in the aligned four-sector class: >5.27x pure-k5 under coefficient norm and >5.49x under fixed uniform-Fisher variance.
- ROBUST-QUARTIC-DRIFT-OPTIMUM-24 certifies the co-phased sign-robust optimum
  `rho_* in [0.71753326806113993,0.71753326806121198]`, >2.222x the pure-k5 robustness radius.
- Full-7D phase optimization repeatedly lands on the same ray, but the rigorous phase inequality for the robustness ratio remains open. Therefore full-7D globality of `rho_*` is **not** promoted.

## B. Dynamic tree inverse lane

**Closed/strengthened:**
- dynamic pole/residue response and terminal-storage intervals under declared minimality/visibility assumptions;
- exact Hankel/Vandermonde conditioning identity: inverse instability separates into static short-edge margin, pole separation and residue visibility;
- explicit two-pole Prony interval certificate;
- maximal positive-storage degree-2 path has constructive Cauer-type uniqueness from the full positive-real two-terminal admittance. Static degree-2 suppression remains an equivalence only when dynamic storage information is absent.

No uniform stable inverse exists without quantitative margins on edges, pole separation and residue visibility.

## C. Hierarchy geometry / partition lane

**Exact algorithms:**
- Gaussian-gauge convolution identity removes the split Boltzmann kernel and permits exact coefficient convolution / FFT implementations.
- exact variance recursion decomposes `C_H` by hierarchy depth.

**8 phases:**
- M=16 peak: `alpha_*≈0.73413`, `C_H≈11.0917`.
- M=32 peak: `alpha_*≈0.74590`, `C_H≈63.1710`.
- M=64 independently validated peak: `alpha_*≈0.82457`, `C_H≈300.610`.
- The apparent early `alpha=beta n` collapse **breaks at M=64**; no stable thermodynamic exponent or transition theorem is accepted.
- Root-level susceptibility fractions at the finite-size peaks: 65.81% -> 77.97% -> 86.62% for M=16,32,64.
- At M=64 the remainder is fully resolved: depth fractions 86.620%, 8.467%, 3.172%, 1.249%, 0.492%.
- Root split distributions are broad/multisector; the peak is not a clean two-state A↔B switch.

**12 phases:**
- all 495 near-equal M=16 compositions have a complete producer;
- all 29 D12 composition orbits were exhaustively analyzed at M=16, and M=32 curvature/selected peaks were computed with exact derivative recursion;
- M=16 pseudocritical points depend strongly on composition (`beta_*` roughly 0.596–0.766 in the full orbit campaign);
- after replication to M=32 the curvature near the common region is much less composition-sensitive, but simple `alpha=beta n` collapse does not hold for the 12-phase 2/1 -> 4/2 sequence;
- root-order distributions at M=32 are broad and multimodal rather than cleanly bimodal;
- an exhaustive M=32 impurity modular rule was obtained for the declared impurity class.

All hierarchy statements remain finite-size/model-conditional. No `D_H=2`, physical spatial geometry or thermodynamic-limit theorem follows.

## D. TWO-MEMORY lane

- Product-model separation is valid but not generic.
- A weak tree->hidden-k2 coupling can imitate the original `k5 nonzero / k3,k4,k6 zero` signature, giving a concrete identifiability counterexample.
- A preregistered phase scan repairs that declared coupling at first order.
- Full first-order hidden-sector tomography uses mixed probes; cubic and quartic channels separate sum/difference information for departure and target perturbations.
- At second order, quartic-only tomography has a no-go: models with identical first-order quartic data can have different O(epsilon^2) offsets.
- Cubic+quartic tomography repairs the explicit second-order obstruction in the declared smooth coupling class and makes the constant k=5 contamination computable from recovered hidden components.

Generic nonlinear coupled identifiability beyond the analyzed orders/classes remains open.

## E. Foundational frontier unchanged

No result here sources physical time, SI length/action, hierarchy scale r, physical storage law, sparse spatial incidence, laboratory evidence, QW-2191, the legacy→strict bridge/role transfer, role-bearing `L_total`, Standard Model, gravity or Theory-of-Everything closure.
