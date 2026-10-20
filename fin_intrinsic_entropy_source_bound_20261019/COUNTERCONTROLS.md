# ID-001 — necessary scope and countercontrols

All of these are *logical control cases*, not physical FIN experiments or replacements for its generator.

1. **Fixed sampling interval vs all resolutions.** For a single δ>0, the discrete-time sampled N2 chain has finite entropy rate hδ. Its two-sided path-space shift is a deterministic invertible transformation with finite entropy hδ, from which the target symbols are read exactly. With a constant-roof suspension of height δ (plus the usual phase coordinate), one obtains a stationary deterministic continuous-time flow of finite time-one entropy proportional to hδ/δ that reproduces that *fixed δ* skeleton. The required H grows as δ→0; there is no contradiction with the theorem. No universal fixed finite H is supplied by this construction.

2. **Discrete-time Bernoulli shifts.** The bilateral shift on i.i.d. q-symbol sequences is invertible deterministic, measure-preserving and has metric entropy at most log q per shift. Deterministic readout of the zeroth coordinate yields an i.i.d. random-looking sequence. Thus 'deterministic source' does **not** imply zero KS entropy or zero stochasticity of symbolic observations.

3. **Infinite-entropy exact path-space source.** The two-sided sample-path shift of the full stationary nontrivial finite-state CTMC (with continuous jump times) defines an invertible, measure-preserving deterministic flow. Its sampled observation yields the exact CTMC at every δ. But h(T1)=∞ for the full clock-resolved CTMC history: any finite-H assumption is violated. This establishes why ID-001 forbids a *source class*, not exact deterministic representation in all senses.

4. **Same stationary p, different dynamics.** Both `B=2(1 p^T−I)` and `B'=c B` for any c≠1 have the same stationary distribution p and generally different jump rates, sampled entropies, and multi-time laws under the **same** model clock. Static p cannot uniquely select activity. (This control does not introduce B' as a FIN law.)

5. **External stochastic readout.** A constant deterministic source has H=0. If one appends fresh independent q12 random draws at each sampling step, its reported Y sequence has entropy rate log12. The implication `h(Y)≤δ H` fails because this is a *randomized*, not deterministic, readout. Supplying a bath, external random clock, or side-channel entropy changes the source class.

6. **Two-time correlations do not imply full stationary joining.** Take independent symmetric signs A,B and define block triples `(A,B,AB)`; concatenate independent blocks and randomize a single three-phase offset to make the resulting sequence stationary. Any two *distinct* time indices are independent unbiased signs, just as in an i.i.d. symmetric-sign process. But `E[Z_0 Z_1 Z_2]=1/3`, whereas it equals zero for i.i.d. signs: only the phase-aligned triple contributes nonzero expectation. Thus all one- and two-time sign moments agree while higher histories differ. Matching a pair correlation or one finite-window TV distance does not furnish the stationary whole-path joining with the needed small per-site mismatch.

7. **Large-error branch of Fano.** The bound `hδ≤δH+F(e)` is valid for all 0≤e≤1. `F_12(e)` is increasing only up to e=11/12, thereafter decreasing. The inverse-error lower bound in the report is therefore stated only for the prescribed branch; no false global inverse or monotonicity is assumed.

8. **QW206 historical 'KS' is a spectral histogram.** The inspected code bins spacings of eigenvalues. It contains neither time-indexed partitions of a flow nor invariant measure nor n→∞ orbit-block entropy; that output cannot be substituted for H in this theorem. Positive genuine KS entropy in an invertible dynamical system does not automatically imply thermodynamic irreversibility.

9. **Same physical time not established.** Replacing the source time δ by an independently tuned δ-specific map can change Hδ and invalidate the common fixed-H conclusion. The flow parameter and N2 heat-bath time here are only *identified as a hypothesis for comparison*.
