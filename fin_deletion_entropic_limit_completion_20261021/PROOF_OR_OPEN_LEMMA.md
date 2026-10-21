# Scope clarification, not a new proof of DEL

Date: 2026-10-21. The complete source derivation remains in its unchanged
[proof map](../fin_deletion_entropic_limit_research_20261020/PROOF_OR_OPEN_LEMMA.md).
All local coefficient, cell, variational, recovery and path arguments are
analytical conditional work by the same author, not independently certified.

The sole extra clarification follows directly from source NC. Put
epsilon_N=E_pi|G_N-dagger-2|->0. Exact stationarity gives

    E sup_(t<=S)|integral_0^t(G_N-dagger-E_pi G_N-dagger)ds|
       <=S E_pi|G_N-dagger-E_pi G_N-dagger|
       <=2S epsilon_N ->0.

Thus the centred IBA condition, as well as the version centred at 2,
HOLDS for this particular comparison. No mixing, irreducibility or
spectral-gap estimate is required for this implication.

It is not a counterexample to the exact-drift KR class: its limiting
phase-conditioned drift is b+partial_z log p_P, not b. All native
drift-matched flow and universal KR IBA obligations remain open. This
clarification is not counted as a second research theorem or a natural
evolution source. No missing lemma is replaced by a PASS label here.
