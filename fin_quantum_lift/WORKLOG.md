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
