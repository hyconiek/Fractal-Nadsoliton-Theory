# New discovery goal after completed ST8651/ST8652

Goal restarted on 9 September 2026. The earlier operational report is
complete and is not counted again as a new discovery.

Selected frontier: test a microscopic many-copy realization of the supplied
normal-ordering source. Tensor replication and a pair-Hamiltonian law are
explicit extra premises; their mathematical consequences are not relabelled
as a physical source derived by FIN.

Candidate exact lift: for X_ij=(|i><j|+|j><i|)/sqrt(2),

    T(rho)=sum X_ij Tr(X_ij rho),
    V=sum X_ij tensor X_ij = (|Omega><Omega|+Swap-2D)/2.

Its Hartree field is T(rho), as required in the previously audited fast-
learning limit. It does not yet reproduce the finite-speed dissipative
controller K_dot, and no controller/noise source has been supplied.

Questions being tested:

1. Does the lift admit exact finite-N quantum dynamics and a controlled
   mean-field limit without the earlier single-state ensemble ambiguity?
2. Can a full-rank strict marginal be a stationary PRODUCT microscopic state?
3. Can the same marginal be stationary with supplied correlations instead?
4. What correlation, entanglement, preparation and symmetry premises are
   actually necessary, rather than merely sufficient in one construction?

Initial code compares the full two-copy unitary, an analytic reduced-state
formula, and an explicit antisymmetric stationary completion. Further proof,
falsification and literature review are required. No PDF generation.

## Final proof audit, 12 September 2026

The main theorem now covers arbitrary nonzero canonical pair couplings on
any finite connected graph, arbitrary added local Hamiltonians, and unequal
full-rank product factors: only the maximally mixed product is stationary.
The exact local commutant and the zero-partial-trace decomposition are the
proof, not numerical stationarity failure. Rank-one product eigenstates and
the n=2 exception are retained as explicit controls.

The same strict marginal has a drifting independent realization and a
stationary antisymmetric correlated completion. The latter is a flat-band
storage construction, not the desired propagator; a common unitary pulse
exposes the distinction. The given sufficient positivity floor cannot be
removed: a pure input makes the affine completion formula nonpositive.

A separate fresh-sample collision processor has an explicit diamond-norm
error bound, an exact one-step uniform-mode leakage formula and a positive
1/N asymptotic leakage coefficient. Program loading is optimized only
within the stated variance-bound architecture. An environment-only example
refutes interpreting variance alone as a universal operational noise cost.
Different programs with the same T can give different finite-copy noise.

Primary mean-field, density-matrix-exponentiation, sample-complexity and
no-programming papers were read. Sample-based simulation is known, and
P512's frequency-ratio theorem is an existing repository result, not a new
discovery. Its application to the new finite rational Hamiltonian class is
kernel sensitive: the canonical legacy cycle has algebraic gap ratios.
Both kernel references separately satisfy the structural correlation tests.

The important new FIN-specific result is the exact correlation requirement
for the declared microscopic source class, together with the separation of
stationary storage from operational propagation. The controlled processor
shows a viable operational route, but still uses a target-encoded program.
No source-independent physical FIN closure or global mathematical priority
is claimed. The final report is
`../FIN_Quantum_Source_Correlation_and_Programming_Report.tex`.

The final scope audit also rejects applying the affine antisymmetric
completion formula at the optimal loading: it becomes nonpositive there.
This is failure of that particular formula, not exclusion of all correlated
completions at the boundary.

The package has 25 scientific tests. `verify.py` records actual replay,
source hashes, static TeX checks and the requirement-by-requirement audit.
No PDF or laboratory record is produced.
