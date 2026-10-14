# FIN KI-001–KI-005 — kinetic invariants and fold response

Date: 2026-10-14. Target directory `fin_kinetic_invariants_campaign_20261014/`. **Source baseline**: user-provided plan, accepted SF intake summary, and hash-preserved predecessor SF/TH archives. No live Git working tree was available for exact HEAD verification; the input source snapshot and mathematical formulas are documented. This is a locally executed *mathematical* campaign, not a repository commit or laboratory experiment.

## Scientific status

**KI-001 — exact conditional identities:** On the 11D occupation tangent, the constrained Gibbs Hessian at the inherited isolated fold has one zero and ten strictly positive directions. With r=Sigma_f Xv and source orientation, the primal coefficients are exactly a_F=a and b_F=b. There is no numerical refit.

**KI-002 — exact kinetic formulas:** HB, MET and g0-matched BAR count generators preserve the same finite-N Gibbs law under symmetric constant-diagonal A7. Their limiting drift admits the logarithmic-mean reversible gradient structure, and at the fold J=-K_f H_F, D_f=2K_f. K_f is positive on the tangent. The soft mobility is m_K=[r^T K_f,T^-1 r]^-1. HB yields m_HB=g_f algebraically; MET and BAR differ while using the same g0 clock.

**KI-003 — prepared compact-stopped analytic limit, conditional on inherited exact fold and its source review:** Using the *actual* microscopic finite-N rates, a 10D stable Lyapunov chart, stopped second energy and a new stopped fourth-moment maximal inequality, the kinetic coordinate x_K and the common thermal readout have the same weak local soft limit under balanced-rounded preparations:

    dZ_K = -m_K [a lambda+(b/2)Z_K^2]ds + sqrt(2m_K)dB_s.

The MET acceptance kink is handled by a locally Lipschitz mobility times a vanishing constrained thermodynamic force, not by falsely invoking smooth C3 rates. The already admitted SF joint-filtration/boundary completion is reused after checking rate/bracket and maximum-fast hypotheses. This is an **author-presented analytic argument**, not a mechanized or independently peer-reviewed stochastic theorem. There is no outward finite-N error rate.

**KI-004 — non-specificity after lawful normalization:** c=(4/b)^(1/3), d_K=c²/(2m_K), Y=Z/c and tau=s/d_K give `dY=(kappa lambda-Y²)d tau+dW`, kappa=-a c/2, for every admitted local lane. The previously rational-certified lambda=0 ±1 eventual upper exit value 0.416855651698382283... therefore is **class-level**, not an A7-specific scalar fingerprint. Under the fixed g0 clock, ratios of mobilities are kinetic fingerprints only.

**KI-005 — scope handoff:** Stop for architect; no new physical source, experimental input, or downstream automatic campaign.

## New numeric diagnostics (binary64 only)

| Kinetic | m_K | m_K/m_HB | d_K=c²/(2m_K) |
|---|---:|---:|---:|
| HB | 3.5156447168395966 | 1.0000000000 | 1.4814354133480332 |
| MET | 0.46077117240094584 | 0.1310630651 | 11.30322575767417 |
| BAR | 0.8694091205061657 | 0.2472972073 | 5.990506036150055 |

Numbers are neither certified intervals nor measurements. The original fold g_f≈3.51564471683959, a≈-0.21257163305224, b≈0.11898259727686 are *reused* established inputs. The independent program rebuilts the finite matrix/fold rather than importing the producer's numerical core; it checked 41 numerical identities/controls. Separately, 396 finite count transition-rate cases passed a producer-side floating detailed-balance check with maximum logarithmic discrepancy about 1.5×10^-14. This is algebraic/numerical corroboration, not an independent stochastic proof.

## Essential proof dependency graph

`accepted kernel/fold/inertia` + `fixed finite-N detailed-balance kinetic rules`
→ `thermodynamic a,b and K_f`
→ `J=-K H_F, one soft and ten stable; m_K, D=2K`
→ `second/fourth stopped moment energy and finite-N LOO control`
→ `joint thermal-readout/kinetic-coordinate path convergence`
→ `same canonical law under kinetic-specific model-time normalization`
→ `negative FIN-specificity for normalized lambda0 ±1 exit outcome`.

The fourth-moment step is indispensable: second-moment tube N^-1/4 does not alone show N^(1/3)*(ell_HB-ell_K)w ->0 uniformly. The proof obtains an epsilon N^-1/3 crossing probability O(epsilon^-4 N^-1/3) for balanced-rounded starts. All constants are existential over frozen finite R/S/lambda; no finite-N certified N0 is asserted.

## Source limitations and no-promotions

The originally admitted exact fold and SF boundary theorem are inherited *conditional* inputs, not reconstructed as an all-operation proof from public source in this run. `source_context/FOLD_SOURCE_BINDING.json` explicitly says its raw certificate bytes were not revalidated during the predecessor run. The previous source helper's BAR coefficient is NOT the present g0-matched BAR convention; the KI plan's 2/12 factor is followed. No changes have been made to AGENTS.md or historical campaigns.

No microscopic kinetic law is derived from equilibrium A7; neither FIN's presence in nature nor physical time, units, space, particles, quantum mechanics, gravity, Standard Model or ToE are established. No global B4 barrier/capacity, metastable waiting time, arbitrary initial preparation, or asymptotic finite-N error rate is inferred.

## Reading order

1. `KI-001/THERMODYNAMIC_MATCHING.md`.
2. `KI-002/REVERSIBILITY_AND_MOBILITY.md`, `KINETIC_RATE_CONTRACTS.json`.
3. `KI-003/KINETIC_CRITICAL_LIMIT.md`, `METROPOLIS_REGULARITY_GATE.md`, `OBSERVABLE_TRANSFER_PROOF.md`, `ESTIMATE_LEDGER.json`.
4. `KI-004/INVARIANT_AND_FINGERPRINT_TABLE.md`, `SPECIFICITY_VERDICT.json`.
5. `KI-005/CLAIM_REGISTER.json` and `NEXT_OPEN_ATOM.md`.
6. `REPLAY.md`, `INDEPENDENT_VERIFICATION.json`, and `MANIFEST.sha256`.
