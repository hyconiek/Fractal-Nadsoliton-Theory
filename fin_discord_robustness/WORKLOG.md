# New goal: microscopic robustness and operational identifiability

Started after the completed separability/necessary-discord report. The first
target was its unresolved exchange-shift value -5/12. This is not another
counting of the prior result.

An exact polynomial reconstruction of V from V+aSwap resolves that value.
Only -1/2 and -3 are genuine stationary-state transfer exceptions at n=12;
explicit positive counterexamples show why they require separate treatment.
Their inverse-block and Hamiltonian norms nevertheless preserve the previous
quantitative non-CQ bound, giving a uniform result for every a,b.

The pure-only class has an exact identity after state swap averaging. It
gives a three-way lower bound involving distance to CQ, state-exchange
asymmetry, and the stationarity residual, without a bound on the arbitrary
antisymmetric Hamiltonian block.

An actual asymmetric CQ equilibrium realizes the average-marginal branch:
its marginals are I/12 and I/12+W/10, not two copies of I/12+W/20. A complete
Luders instrument prepares it from one supplied amplified program. Passive
swap-covariant records agree with its symmetrized canonical description,
but labelled marginal measurements and coherent controlled-U distinguish
them. The latter counterexample prevents an overbroad operational claim.

Two universal stationary projector channels have identical complete local
channels. Their mixture is separable exactly at weight 1/2 and NPT otherwise,
for every input state. The balanced channel uses one program as a global
CPTP map; its negative marginal Choi partial transpose rules out universal
LOCC broadcasting from one unknown input held at a fixed party without an
additional quantum resource. This does not prohibit local preparation of a
fixed known separable target from its classical description.

Primary sources on controlled unknown operations, reference-frame restrictions
and quantum cloning were inspected. Their general mechanisms are not claimed
as new discoveries. No physical source, calibration, apparatus or ToE closure
is exported. Output policy remains source-only TeX for the final report; no PDF.

## Final verification scope

The package has 20 scientific tests and one replayable result set. Exact
polynomial/band identities, symbolic marginal identities and partial-
transpose minors are the proof layer. Numerical checks separately test
cross-band exceptions, passive records, coherent-control discrimination,
and the stationary channel family. A zero-program control shows that
entanglement alone does not create a nonzero local kernel. A separate
legacy program verifies the kernel-independent channel identities without
transferring strict numerical discord bounds or physical roles.

The single final report is
`../FIN_Discord_Robustness_and_Operational_Identifiability.tex`.
`verify.py` records execution, source hashes and the completion audit. The
remaining physical obligation is a sourced joint preparation law and
operational access model, not a choice inferred from local tomography alone.
