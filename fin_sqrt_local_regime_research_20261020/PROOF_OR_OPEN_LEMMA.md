# Terminal proof and still-open dynamic gate

Date: 2026-10-20. **CONDITIONAL_LOCAL_WINDOWS_PROVED_BY_COMPARISON.**

The proof is split into three fully specified steps:

1. [Exact rates and conductances](RATE_AND_FORM_COMPARISON.md): Gibbs balance,
   bounded rates and a_* E_HB<=E_SQRT<=a^* E_HB on the actual MT domain.
2. [SQRT local window](LOCAL_WINDOW_THEOREM.md): invariant-weight domination,
   departure compensator, actual iid density and pre-exit coupling give
   unpostselected local equilibrium and one-sided exponential residence.
3. [Finite-source window](FINITE_SOURCE_DIRECT_THEOREM.md): c_N>=2(1+RN)
   pays rate comparison and discrete laziness. Direct stationary-weight
   domination gives the same property for an actual finite entropy source,
   without a continuous-path approximation or exchanging limits.

These proofs inherit actual MT's static geometry, all-type connectivity,
HB gap and moment results. They are a new application to an already
derived law and a direct finite-source corollary, not new proofs of those
old hypotheses. See INPUTS.sha256 for the exact files used.

## First remaining process question, not a missing part of this theorem

A stationary critical path theorem for SQRT needs its own parity drift,
fast-direction selection, bracket, tightness and joint-filtration argument.
The formal expected coefficient is q times the HB critical coefficient,
because the new zero-field per-copy rate is 1, not 1/q. Neither Dirichlet
comparison nor the static quartic limit proves it. Exact native statewise
parity drift is different; the KR exact-drift theorem is not applicable.

For a finite source at c_N=O(N), the factors r_x r_y in its true rates need
not tend to 1. Its critical coefficient cannot simply be called q. A
separate scaling with c_N/N->infinity must pay the SL horizon/error budget,
or a native proof must retain and control those factors. No such critical
trajectory theorem is claimed in this packet.

The next finite instruction is [NEXT_CRITICAL_ORDER.md](NEXT_CRITICAL_ORDER.md).
Physical sourcing of the entropy principle, supplied Gibbs/A7/g/carrier,
preparation and comparison-to-physical clock remains open in either case.
