# FIN BA-001–BA-004: exact parity-cut rigidity; universal bracket averaging still open

**Date:** 2026-10-16. **Origin:** user-provided post-KR BA order and KR acceptance summary. **Status:** `EXACT_PARITY_COMPRESSION_RIGIDITY_WITH_RESTRICTED_VARIATIONAL_NEGATIVE_CONTROL; Q12_UNCONDITIONAL_IBA_OPEN`.

This handoff contains **only new BA work**. The primary q12/A7 heat-bath, conditional MT Gibbs state, KR-002 critical theorem under NC/IBA, all predecessors and `AGENTS.md` are unchanged. This is intrinsic model mathematics, not a newly sourced physical kinetic law.

## Precise new positive theorem

For every reversible diagnostic count generator on the actual reflected MT domain, using only single-copy replacements, sharing exactly the native reflected scalar parity drift `d(x)` at **every** count state and stationary Gibbs `pi_D,N`, the total equilibrium conductance crossing each adjacent parity level is *exactly determined* by `pi_D,N` and `d`. Specifically, with `delta=2/N`,

`K_l= (1/delta) sum_(x:m(x)<=a_l) pi_D,N(x)d(x)`.

Hence the parity-compressed instantaneous generator, every bilinear Dirichlet form on parity functions and the **entire parity-conditional mean of the microscopic bracket** are exactly invariant under arbitrary admissible microscopic activity redistribution:

`E_pi[G_tilde|m=a_l]= E_pi[G_native|m=a_l]`, `G=N Gamma(m,m)`.

Native UF continuity plus MT equilibrium imply `||E_pi[G_tilde-E_pi G_tilde |m]||_L1(pi)->0` and, since `0<=G<=4C`, also L2 convergence. This is stronger than the previous smooth static test from KR and needs neither NC nor IBA. It fixes conditional *means*, not within-fibre fluctuations or their history.

## Exact component obstruction

The positive-rate graph need not be connected even though the original MT domain is connected. With `C` a diagnostic communicating component, `h=P_kerL(G-E G)` obeys

`h|C=-2N E_pi[m d|C]+2N E_pi[m d]`.

The formula is exact, and the accepted hypotheses do not give an all-N bound on `||h||_pi`. If all permitted edges are strictly positive, `h=0`, but there is no uniform lower activity or Poisson bound from that fact.

## BA-002 Poisson/Dirichlet attempt

The finite reversible forward/backward martingale decomposition yields exactly

`E_pi sup_(t<=S)|int_0^t r(X_(sqrt(N)s))ds|² <= 2S²||h||²_pi+40S R/sqrt(N)`,

where `r=G-E G`, `-Lu=r-h` on components, `R=<r-h,u>`. Thus `||h||->0` and `R=o(sqrt(N))` imply IBA without NC. These estimates are **not** established uniformly over the actual q12 diagnostic conductance polytope.

New six-state rational negative control: all equilibrium parity cut traffic and all statewise scalar drifts stay unchanged under an explicit wheel stress; nevertheless microscopic Gamma changes and the limiting stressed chain has nonzero invariant-component noise means. Along a strictly positive approach with intercomponent conductance `tau`, its exact Poisson cost diverges as `R(tau) ~ (1064/126075)/tau`. This disproves a generic shortcut from drift/Dirichlet projection to uniform mixing. It is **not** an admissible q12 critical counterfamily, and it does **not** refute universal q12 IBA.

## Terminal BA-003 decision

No all-N q12/A7 counterfamily with a different limiting bracket was constructed; no class-specific variational upper bound or direct time averaging proof without NC was completed. The missing lemma is precisely the conditional fibre/invariant-component averaging estimate on the given MT count graph under the exact statewise drift flux constraints. Preserve `OPEN_Q12_UNCONDITIONAL_IBA`.

## Source provenance

The available KR source work, MT and UF actual archived files, and a historical AGENTS snapshot were hashed. The *current* post-KR review folder, complete activity-covariance research files, and current checkout were not available. Consequently no declaration of independent verification of those original bytes is made. The supplied intake and proof criterion were used as user-provided mathematical premises; this handoff is suitable for independent audit with current repository inputs. See `SOURCE_PROVENANCE.json` and `INPUTS.sha256`.

## Independent diagnostics and replay

`python3 verify_ba.py` runs **24** exact rational checks with SymPy (no Monte Carlo, no large N sweep). Six-state fixture stresses are checks of algebraic constraints, **not** a substitute for the full q12 analysis. `python3 verify_package.py` validates all manifest entries and bounded JSON. No previous package was edited or reserialized.

## One next atom and stop

Either pay a q12 class-specific time-integrated bracket estimate without NC, or construct a full sequence of actual reflected q12 conductances obeying every frozen hypothesis and a genuinely different limiting bracket/history. A huge Poisson dual norm by itself is insufficient to declare failure of IBA.

**Physical limitations:** no natural source of A7 or activity, physical clock/units, spatial incidence, locality, conserved physical energy, apparatus/EV, ontology/selector/role transfer, QM/GR/SM/L_total or ToE. The independent AS/KB physical activity-provider gate remains OPEN.
