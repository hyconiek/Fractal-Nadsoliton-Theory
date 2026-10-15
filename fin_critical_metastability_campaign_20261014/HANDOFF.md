# FIN MT-001–MT-005 — Critical local trapping, conditional Gibbs and finite-N mixing

Date: 2026-10-14. Execution: `fin_critical_metastability_campaign_20261014/`.

**SCIENTIFIC VERDICT: `LOCAL_CRITICAL_WINDOW_PROOF_SUBMITTED_FOR_INDEPENDENT_INTAKE`.**

The proposed asymptotic window is mathematically stronger than the incoming UF finite-horizon result, but its proof is an **authorial analytical submission**, not intake acceptance, proof-assistant certification or independent peer review. The key nontrivial new lemma is the *energy-controlled root-path construction* of MT-001/MT-004. If an independent reviewer finds a defect there, retract the polynomial gap/window claim while preserving independently valid trapping and conditional Gibbs statements.

## Scope and theorem statement

In the same q=12 exact-decimal A7 generalized Curie–Weiss FIN parent, with unit-per-copy **actual leave-one-out** heat-bath dynamics and iid-uniform labels at t=0, let `g_N=g_c+lambda/sqrt(N)` with `g_c=12/Lambda6` exact, fixed bounded lambda and N divisible by12. The two already frozen MET/BAR controls admit the same reversible-weight argument at their specified g0 rate normalization; they are diagnostic, not replacement FIN axioms.

Define ε=1/100, fast radius δ=1/10000, `h=10^-10`, `b=5·10^-11`; let D_N be the count-graph component containing uniform counts of the strict local cylinder `|m|<ε, ||w||<δ` and the sublevel `F_gN(p)≤h`. Under the accepted exact spectral gap at `g_c` and the local analytic lemmas in MT-001, for large N∈12ℕ:

* **MT-001:** D_N is connected and contains iid-uniform starts with probability 1-o(1); all outgoing departure types have `F_gN≥b`.
* **MT-002:** `P_iid(tau_D≤T)≤mu0(D_N^c)+2NT(N+1)^12 exp(−Nb)` for all T≥0. In particular it is o(1) at T=exp(cN), any 0<c<b.
* **MT-003:** Exact conditional count Gibbs `pi_D,N` has joint `(N^(1/4)m, sqrt(N) E_fast^T w)` limit equal to the quartic density in Z times a ten-dimensional Gaussian with precision `E_fast^T(12I−g_c A)E_fast`. All fixed polynomial moments converge.
* **MT-004:** A direct energy-monotone/central-log-corridor **canonical-path** construction gives `gap(reflected D_N)≥c*(N+1)^(-K)` for some fixed K. The complementary Rayleigh bound `gap≤C N^(-1/2)` shows no constant critical gap is claimed. A valid but extremely loose `t_loc=(N+1)^(K+2)` therefore obeys, uniformly in `t_loc≤t≤exp(cN)`, `TV(Law_iid[p_N(t)],pi_D,N)→0`. Both unmatched timescales and the asymptotic `N0` are unquantified.
* **MT-005:** Preserve the full scope, arithmetic statuses, negative controls and stop for the architect.

The reflected chain is an analytical comparison process only; **the parent in the TV theorem is the unaltered, unconditioned chain**, and outside/exit events remain in its TV error. Neither full-Gibbs quartic equilibrium, global mixing/escape, parity-bit lifetime, scalar physical seconds, source-selection of A7 nor independent physical observation follows.

## Files and logical proof map

Read `MT-001_COERCIVITY_CONNECTIVITY.md` first. Its local-gradient/discrete-descent/central-corridor estimates are exactly what license the polynomial gap and cannot be omitted. Then read `MT-002_IID_TRAPPING.md`, `MT-003_CONDITIONAL_GIBBS.md`, `MT-004_REFLECTED_MIXING_WINDOW.md`, and `NEXT_OPEN_ATOM.md`. The five `MT-00X/` directories hold separate reports, structured results, replay and input hashes. `mt_analysis.py` and `independent_verify.py` provide *tiny numerical/rational* checks, not stochastic theorem proofs. Original inputs and UF predecessor are in `source_context/` (unchanged copies). `RESEARCH_SUMMARY_PL.md` presents the Polish summary.

## Reproducibility and custody

Input hashes in `PROVENANCE.json`; generated output hashes in `MANIFEST.sha256`. Neither a producer Git checkout nor all repository review artifacts were available for inspection; the uploaded UF intake is a supplied scope/accepted-premise record. Spectral gap `gamma>0.0186` and exact g_c are **inherited input premises**, not newly certified by the binary64 program. Rational domain comparisons are exact under that premise. No finite-N numerical N0, gap prefactor, polynomial exponent optimised for hardware or achieved experiment is claimed.

**Next:** independent architecture/intake audit of the discrete root-path lemma and lattice dominated-convergence argument, followed by one architect decision. Do not automatically activate B4, FIN-M0, task335, new kernel/selector, physical graph/units or empirical lane.
