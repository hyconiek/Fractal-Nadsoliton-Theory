# RELEASE 8 STRICT TEXTBOOK EDITION (EN + PL)

**Version:** 8.0.0 (Current-state projection — O3 strict gate closure packetized)  
**Date:** 2026-05-12  
**Branch:** `main`  
**Predecessor:** Release 7 Strict Textbook Edition

## Status discipline (no false pass)

- This Release 8 document describes the **current repo state only**.
- It records that the O3 pipeline (`P1313..P1339`) reaches strict closure **inside the declared pipeline rigor** (`P1338` + `P1339`), with explicit machine reports.
- It does **not** claim global ToE closure in the strongest sense.
- It does **not** erase the strict/non-strict distinction in historical packets; non-strict packets remain explicitly labeled.
- It does **not** silently transfer legacy kernel roles into strict kernel claims beyond what is explicitly exported.

---

## 1. Release‑8 target (EN): “Current-state strict O3 closure snapshot”

Release 8 is a state report whose target is:

1. freeze the current strict O3 gate result for `QW-2191`,
2. show exactly which obligations were passed and how,
3. preserve hard limits and scope boundaries,
4. provide a clean handoff to post-closure regression/replication governance.

This release is **not** a restart of theory-building; it is a disciplined snapshot of the executed `P1313..P1339` chain.

---

## 2. Current strict state entering Release 8 (EN)

On the current repo state:

1. O3 plan and step chain are fully packetized (`P1313..P1339`) with explicit pass/fail discipline.
2. Candidate-class and normalization layer are frozen (`P1314`, `P1315`).
3. Pairwise and sweep/replay stages were executed with artifacts (`P1316`, `P1317`, `P1318`).
4. Residual-slot neutrality initially failed (`P1319`), forcing further construction.
5. Non-strict closure lane under explicit premise was separated and labeled (`P1320`).
6. Candidate evolution proceeded (`v2` falsified, `v3` degenerate, `v4` robust under revised gate) (`P1323..P1326`).
7. Formal export lane progressed from `1/2` to `2/2` via final global-L2 attempt (`P1328`, `P1338`).
8. Full checker refresh reached `5/5` obligations (`P1337` refreshed after `P1338`).
9. Independent post-closure replication packet passed (`P1339`).

Current pipeline-level status:

- `o3_strict_ready = true`
- `qw2191_strict_status = CLOSED`

within the declared Release-8 strict pipeline semantics.

---

## 3. What Release‑8 still does not “false‑PASS” (EN)

Release 8 still does **not** claim, unless separately exported:

1. strongest-form global ToE closure beyond this O3 gate scope,
2. retrospective reinterpretation of previously labeled non-strict packets as strict,
3. silent replacement of theorem-level export by numeric success alone,
4. immunity from future independent falsification.

---

## 4. Release‑8 deliverables (EN)

1. This Release-8 textbook snapshot (`RELEASE_8_STRICT_TEXTBOOK_EDITION_EN_PL.md`).
2. Packet chain `P1313..P1339` as the current strict O3 trail.
3. Generated machine reports under `fundamental_action_reconstruction/generated/` documenting each step.
4. Post-closure next-step target: regression/replication guardrail packet (`P1340`).

---

## 1. Cel Release‑8 (PL): „Migawka aktualnego stanu strict po O3”

Release 8 ma jeden cel: **uczciwie opisać obecny stan repo** po wykonaniu łańcucha `P1313..P1339`, bez dopisywania nowych roszczeń ponad to, co jest wyeksportowane.

---

## 2. Aktualny stan strict na wejściu do Release 8 (PL)

Na obecnym stanie repo:

1. Plan O3 i pełna ścieżka kroków są jawnie zapakowane (`P1313..P1339`).
2. Warstwa klas kandydatów i normalizacji jest zamrożona (`P1314`, `P1315`).
3. Pairwise + sweep/replay zostały wykonane i udokumentowane (`P1316`, `P1317`, `P1318`).
4. Pierwotna próba neutralności residual-slot nie przeszła (`P1319`), co wymusiło dalszą konstrukcję.
5. Tor non-strict pod jawną premise został oddzielony i oznaczony (`P1320`).
6. Kandydaci ewoluowali: `v2` obalony, `v3` zdegenerowany, `v4` przeszedł poprawioną bramkę (`P1323..P1326`).
7. Tor formalnych eksportów przeszedł z `1/2` do `2/2` przez finalną próbę globalnego L2 (`P1328`, `P1338`).
8. Pełny checker po odświeżeniu pokazuje `5/5` obowiązków (`P1337` po `P1338`).
9. Niezależna replikacja post-closure przeszła (`P1339`).

Bieżący status w semantyce tego pipeline:

- `o3_strict_ready = true`
- `qw2191_strict_status = CLOSED`

---

## 3. Czego Release‑8 nadal nie może „false‑PASS” (PL)

Release 8 nadal **nie** rości bez osobnego eksportu:

1. najsilniejszego globalnego domknięcia ToE poza zakresem tej bramki O3,
2. cofnięcia etykiet non-strict w starszych pakietach,
3. zastąpienia eksportu theorem-level samym wynikiem numerycznym,
4. gwarancji, że przyszły niezależny audyt niczego nie podważy.

---

## 4. Deliverables Release‑8 (PL)

1. Ten dokument Release 8 (`RELEASE_8_STRICT_TEXTBOOK_EDITION_EN_PL.md`).
2. Łańcuch pakietów `P1313..P1339` jako aktualny strict O3 trail.
3. Raporty maszynowe w `fundamental_action_reconstruction/generated/` dla każdego kroku.
4. Następny krok po domknięciu: pakiet guardrail regresji/replikacji (`P1340`).

---

## 5. QW-2191 closure and physics-consistency audit (EN+PL)

### 5.1 Strict statement boundary

Under the current exported chain (`P1313..P1339`) the repository exports a **pipeline-strict** closure status for the O3 gate:

- `qw2191_strict_status = CLOSED` (inside declared checker semantics),
- with explicit reminder that this is not equivalent to strongest global ToE closure.

This is physically consistent with current guardrails **only if** we keep the closure in the selector/uniqueness-gate scope and do not overclaim host-level SM/GR derivation from O3 alone.

### 5.2 Core equations used in Release 8

Strict operational kernel:

\[
K_{\text{strict}}(d)=\frac{\cos(\omega d+\phi)}{1+\beta d^{\eta}}
\]

Candidate-v4 score map:

\[
s(\varphi,a)=(\varphi-0.40)+0.9(a-0.58)-0.6(\varphi-0.40)^2
\]

Sign-selection rule:

\[
\operatorname{sign}(s)=
\begin{cases}
+1, & s\ge 0,\\
-1, & s<0.
\end{cases}
\]

### 5.3 Derivation history (compressed, release-ready)

1. **Kernel split made explicit** (`K1`, `K2`, `F2`): strict kernel is an operational working kernel; no silent transfer of legacy role-claims.
2. **Frontier classification** (`F3`): continue on kernel-split-robust lane, keeping QW-2191 as a real strict-core obstruction unless exported evidence is produced.
3. **O3 construction chain** (`P1313..P1339`): candidate classes, normalization, pairwise/sweep/replay, residual-slot handling, candidate evolution (`v2->v3->v4`), theorem-level export obligations, global-L2 completion, independent replication.
4. **Exported endpoint**: checker obligations are satisfied in the declared pipeline semantics and the O3 gate status is set to strict closed.

### 5.4 Physical reading discipline

- O3 closure contributes to **selector uniqueness discipline** and stabilizes one admissible branch policy in the strict operational lane.
- O3 closure **does not by itself** derive full SM+GR dynamics; it is one gate in a larger architecture.
- Any claim of full host-identification must still be made via explicit theorem-level exports outside this gate.

### 5.5 PL — krótki werdykt fizyczny

Domknięcie `QW-2191` w Release 8 jest fizycznie spójne **w granicach semantyki pipeline O3**: usuwa niejednoznaczność selektora, ale nie zastępuje pełnej, niezależnej identyfikacji wszystkich sektorów SM/GR.

---

## 6. Candidate strict-core-to-SM effective Lagrangian scaffold (release note)

For release documentation completeness, we record the strict-core candidate density (nadsoliton carrier + 12-octave sector):

\[
\mathcal{L}_{\text{core}} = \frac{1}{2}\partial_\mu\phi\,\partial^\mu\phi
+\sum_{i=0}^{11}\left[\frac{1}{2}\partial_\mu\psi_i\,\partial^\mu\psi_i - V(\psi_i)\right]
+\frac{1}{2}\sum_{i\neq j}K_{ij}\,\psi_i\psi_j,
\]

with
\[
K_{ij}=K_{\text{strict}}(d_{ij})=\frac{\cos(\omega d_{ij}+\phi)}{1+\beta d_{ij}^{\eta}}.
\]

A gauge-covariant extension candidate toward SM-shaped sectors may be written schematically as:

\[
\mathcal{L}_{\text{eff}}=\mathcal{L}_{\text{gauge}}+\mathcal{L}_{\text{fermion}}+\mathcal{L}_{\text{Higgs}}+\mathcal{L}_{\text{Yukawa}}+\Delta\mathcal{L}_{\text{nadsoliton}},
\]

where:

- \(\mathcal{L}_{\text{gauge}}=-\tfrac14 G^A_{\mu\nu}G^{A\mu\nu}-\tfrac14 W^I_{\mu\nu}W^{I\mu\nu}-\tfrac14 B_{\mu\nu}B^{\mu\nu}\),
- \(\mathcal{L}_{\text{fermion}}=\sum_f \bar\Psi_f i\gamma^\mu D_\mu\Psi_f\),
- \(\mathcal{L}_{\text{Higgs}}=(D_\mu H)^\dagger(D^\mu H)-V(H)\),
- \(\mathcal{L}_{\text{Yukawa}}=-y_u\bar Q\tilde H u_R-y_d\bar Q H d_R-y_e\bar L H e_R+\text{h.c.}\),
- \(\Delta\mathcal{L}_{\text{nadsoliton}}\) encodes strict-core residual couplings to emergent sectors.

This section is explicitly a **structured candidate scaffold**, not a claim that full SM matching is already discharged by O3 closure alone.
