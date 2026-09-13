# FIN — pełny handoff badań wykonanych w tej rozmowie
## Od badań poprzedzających `FIN_Discord_Robustness_and_Operational_Identifiability` do bieżącego frontu rank‑7 / Morse-index

**Data handoffu:** 2026-09-13  
**Cel:** przekazanie drugiemu agentowi AI kompletnego pakietu nowych badań wykonanych w tej rozmowie, łącznie z badaniami **sprzed** raportu `FIN_Discord_Robustness_and_Operational_Identifiability`, wynikami samego punktu kontrolnego oraz wszystkimi badaniami późniejszymi.  
**Zakres:** matematyka FIN / strict finite operator / conditional active-gain model / rank‑7 mediator / dual, fazy, Morse, curvature ceiling, Ising/CRT, parity decomposition i rank-one resolvent.  
**Nie jest to nowy oficjalny release repo.** To handoff roboczy z dokładnym ledgerem statusów.

---

# 0. Epistemic legend — obowiązkowo zachować

W dalszej części używane są statusy:

- **[STRICT / REPO THEOREM]** — wynik już obecny/certyfikowany w repo lub wynik dokładny wynikający z frozen strict operator.
- **[ANALYTIC / EXACT FOR SUPPLIED DECIMALS]** — algebraiczny wynik dla dostarczonych wartości liczbowych spektrum; do pełnego theorem-grade upstream nadal warto wykonać outward-rounded interval lift samych `lambda_k`.
- **[COMPUTER-ASSISTED / NUMERICAL CERTIFICATE]** — skończony cover/Bernstein/Krawczyk lub bardzo mocna skończona certyfikacja, ale trzeba zwrócić uwagę, czy współczynniki wejściowe są intervalowe czy tylko dziesiętne.
- **[STRONG NUMERICAL EVIDENCE]** — multistart / differential evolution / Sobol / continuation bez pełnego globalnego certyfikatu.
- **[CONDITIONAL MODEL]** — wynik dotyczy dostarczonego aktywnego funkcjonału lub innej jawnie dostarczonej konstrukcji; nie jest źródłem fizycznym FIN.
- **[REFUTED]** — hipoteza sprawdzona i obalona; nie wracać do niej bez nowej idei.

Najważniejszy guardrail: active gain `g`, jego znak i jego wielkość **nie są strict-derived**. Repo ST293 dowodzi no-go dla pasywnej klasy `C=f(A)>=0`: nie może ona wygenerować ujemnego aktywnego członu informacji. Active/state-dependent/pumped/nonnormal/nonequilibrium nadal pozostają otwarte. Pełny strict operator przy `g=4` ma theorem globalnych 12 minimów w jednym orbicie D12 (ST389/ST390), ale ten theorem **nie przenosi się automatycznie na rank‑7** (ST391/ST410).

---

# 1. Frozen strict data używane we wszystkich nowych obliczeniach

Strict weighted circulant `K12`, shell weights:

\[
w_1=0.4699856726450201,\quad
w_2=0.1920435516901028,\quad
w_3=0.09142861427792495,
\]

\[
w_4=0.0470291687456504,\quad
w_5=0.02413122336363006,\quad
w_6=0.011070817321442113.
\]

`A = diag(row sums) - W`.

Dodatnie eigenvalues Laplacianu / Fourier sectors:

\[
\lambda_1=0.7541211542070795,
\]
\[
\lambda_2=1.5770495144276093,
\]
\[
\lambda_3=1.9614068619764449,
\]
\[
\lambda_4=2.199568849333209,
\]
\[
\lambda_5=2.2986062720790956,
\]
\[
\lambda_6=2.3421820411463004.
\]

Multiplicities: `(2,2,2,2,2,1)`. Mean-zero rank `11`.

Uniform state:

\[
u=\frac1{12}{\bf 1}.
\]

Conditional active-gain functional:

\[
V_g(p)=D(p\|u)-\frac g2(p-u)^TA(p-u).
\]

Uniform local spinodal dla pełnego `A`:

\[
g_{\rm spin}=\frac{12}{\lambda_6}
=5.1234275513986\ldots
\]

Repo: active sign i `g` supplied, nie wyprowadzone.

---

# PART I — BADANIA WYKONANE W TEJ ROZMOWIE **PRZED** RAPORTEM DISCORD

# 2. Relation–relation / edge-current structure

## 2.1 Exact Hodge decomposition of edge-current noise

Dla kompletnego weighted graph `K12`:

- liczba krawędzi: `66`,
- wymiar gradient/image incidence: `11`,
- wymiar cycle space: `55`.

Kanoniczna cykliczna kowariancja prądów:

\[
\Sigma_{\rm cyc}
=
G-GDA^+D^TG
=
G^{1/2}P_{\rm cycle}G^{1/2}.
\]

Wyniki:

- rank `55`,
- `28` różnych dodatnich eigenvalues,
- \(\lambda_{\min}^+\approx 0.01250017\),
- \(\lambda_{\max}\approx0.46998567\),
- condition number \(\approx37.5983\),
- D12 covariance residual \(\approx5.5\times10^{-15}\).

Dokładny rozkład całego szumu krawędziowego:

- rank części zmieniającej `p`: `11`,
- rank części czysto cyklicznej: `55`,
- `||D^T C_cycle|| ~ 8.4e-15`,
- state-changing contraction reprodukuje `A` do ~`2.1e-14`.

**Interpretacja:** divergence-free cycle currents są dużą relacyjną przestrzenią pamięci/ukrytego ruchu, ale liniowo nie zmieniają `p`. To nie jest dodatkowy „byt pod nadsolitonem”; jest to struktura relacji wewnątrz całej konfiguracji.

## 2.2 Crossing geometry is not intrinsic

Testowano pomysł, czy przecięcia linii relacji w konkretnym rysunku mogą stanowić dodatkowe „interakcje relacja–relacja”. Po permutacji/innym osadzeniu tych samych krawędzi liczba geometrycznych crossings się zmienia.

**Wniosek:** crossings konkretnego planarnego rysunku nie są intrinsic FIN object, jeśli nie zostanie dostarczona dodatkowa geometria/embedding law.

## 2.3 Passive independent-walker model

Dla `N` niezależnych walkerów na strict graph:

\[
\dot p=-Ap,
\]

a lokalna noise covariance:

\[
\operatorname{Cov}(dp)=\frac{A}{6N}\,dt.
\]

Modalnie:

\[
d q_k = -\lambda_k q_k\,dt + \sqrt{\lambda_k/(6N)}\,dW_k.
\]

Stationary variance:

\[
\operatorname{Var}(q_k)
=
\frac{\lambda_k/(6N)}{2\lambda_k}
=
\frac1{12N}.
\]

**Wniosek:** strict diffusion seeduje wszystkie mody, ale damping rośnie dokładnie z tym samym \(\lambda_k\); noise + damping równoważą się. Passive walkers **nie generują active gain ani localization**.

To jest ważny no-go dla intuicji „sam szum wybierze fazę”.

---

# 3. Passive strict memory vs active gain

Even/odd Schur reduction strict operator:

- hidden dimension `6`,
- minimal visible realization dimension `5`,
- tylko `3` visible pole groups,
- `Tr Sigma(0) ≈ 2.055115`,
- hidden eigenvalue range ~`1.171091 ... 1.961407`.

Schur self-energy jest typu Stieltjes / passive memory. Efektywny operator zachowuje nieujemność:

- `min eig A_eff ~ 1.5e-16`,
- zero negative eigenvalues,
- one zero mode,
- remaining spectrum positive.

**Wniosek:** strict Schur memory może obniżać stiffness, tworzyć retardation/pole structure i ukrytą pamięć, ale nie odwraca znaku do active negative-information gain.

To zgadza się z repo ST293. Nie wolno używać pasywnej pamięci jako „ukrytego źródła g”.

## 3.1 Mean-field multi-copy realization

Zbadano warunkową realizację typu mean-field/Gibbs, w której `V_g` może pojawić się z parametrem:

\[
g=\beta J.
\]

To pokazuje, że funkcjonał ma sens jako free-energy w pewnej klasie modeli. Ale `beta` i `J` są dodatkowymi zasobami/model parameters.

**Status:** [CONDITIONAL MODEL], nie strict provenance.

---

# 4. Spectral rank ladder — dlaczego rank 7 jest specjalny energetycznie, ale nie informacyjnie

Rank‑7 mediator = retained top Fourier blocks `k=3,4,5,6` (+ conjugates), total rank:

\[
2+2+2+1=7.
\]

Equivalently low sectors `k=1,2` są usunięte.

Real-space kernel `A7` jest signed high-pass; przykładowy pierwszy wiersz:

- `d0`: `+1.27177883399365`
- `d1`: `-0.71025447835975`
- `d2`: `-0.12346618833839`
- `d3`: `+0.17141297146001`
- `d4`: `+0.14723505364057`
- `d5`: `-0.04670400338684`
- `d6`: `-0.14822554402486`

To tłumaczy, dlaczego positive-weight cap/rearrangement methods dla pełnego `A` nie transferują.

## 4.1 Exact spectral-budget threshold at rank 7

Dla D12-covariant top-block `B <= A`:

\[
V_g(e_j)-V_g(u)
=
\log 12 - \frac{g}{24}\operatorname{Tr}B.
\]

Przy `g=4` vertex może pokonać uniform tylko gdy:

\[
\operatorname{Tr}B>6\log 12
\approx14.9094398987.
\]

Top rank-6 budget:

\[
S_6=2(\lambda_5+\lambda_4+\lambda_3)
\approx12.919164
\]

— niewystarczający.

Top rank-7:

\[
S_7=S_6+\lambda_6
\approx15.261346
\]

— wystarczający.

Margin:

\[
S_7-6\log12\approx0.351906.
\]

**Exact conclusion:** rank 7 jest minimalnym D12-covariant **energetic rank** w tej top-block family, który przy `g=4` pozwala pure vertex pokonać uniform.

Nie oznacza to, że rank 7 jest minimalnym wymiarem „informacji” ani minimalnym wymiarem rozróżniania 12 etykiet.

## 4.2 Vertex signatures: distinguishability arrives earlier

Cumulative top blocks:

- rank 1 (`k=6`) daje tylko `2` signatures;
- rank 3 (dodanie pair `k=5`) daje już `12` distinct vertex signatures;
- rank 5,7,9,11 też 12.

Minimal signature distance:
- rank3 ~`0.408248`,
- rank5 ~`0.816497`,
- rank7 `1`,
- rank11 `1.414214`.

**Wniosek:** rank 3 wystarcza do identyfikacji 12 labels w feature space. Rank 7 jest progiem energetycznej localization przy `g=4`, nie magiczną liczbą informacyjną.

## 4.3 Nonlinear Fourier closure

Aktywne Fourier modes generują nieaktywne przez exponential/softmax mixing.

Cumulative rank:
- rank1: tylko modes `0,6`;
- rank3: wszystkie 12 modes generowalne do order 3;
- rank5: wszystkie do order 2;
- rank7: wszystkie do order 2.

To jest ważne dla późniejszego exact dual.

---

# 5. Numerical transition ladder vs active rank

Strong numerical multistart thresholds:

| active rank | numerical first competitor `g_c` | `p_max` |
|---:|---:|---:|
| 1 | 5.123428 | ~0.08334 |
| 3 | 5.123428 | ~0.08335 |
| 5 | 4.821627 | 0.667433 |
| 7 | 3.7183448981 | 0.836365 |
| 9 | 3.121905 | 0.884995 |
| 11 | 2.902496 | 0.904841 |

Pure-vertex threshold dla rank7 ~`3.907765`, więc first localized competitor pojawia się przed pure vertex.

Pełny strict rank11 ma repo theorem/certificates; rank7 numerical/global proof nadal open.

---

# 6. Exact 11D -> 7D dual and max-entropy completion

Dla rank7:

\[
A_7=XX^T,
\qquad X\in\mathbb R^{12\times 7}.
\]

Dual:

\[
\Phi_g(\theta)
=
\frac{\|\theta\|^2}{2g}
-
\log\left[\frac1{12}\sum_{i=1}^{12}e^{\theta\cdot X_i}\right].
\]

Stationarity:

\[
h=X\theta = gA_7(p-u),
\]

\[
p_i=
\frac{e^{h_i}}{\sum_j e^{h_j}}.
\]

**Exact conceptual result:** aktywne pole `h` jest 7D i ma zero components w omitted Fourier sectors `k=1,2`, ale pełny stan `p` może mieć duże low-mode content wskutek nonlinear exponential mixing.

Przy rank7 coexistence około `32.9075%` power `q=p-u` leży w inactive modes.

Jeśli sztucznie wymusić:

\[
p-u \in {\rm Range}(A_7),
\]

transition threshold rośnie z:

\[
g_c^{(7)}\approx3.7183448981
\]

do:

\[
g_c^{(7,linear)}\approx4.720214.
\]

Reduction:

\[
\Delta g\approx1.001869
\]

(~21.2%).

**Interpretacja:** w conditional model brakujące linear coordinates są uzupełniane przez max-entropy completion. Nie są „dodatkowym aktywnym polem”; powstają jako nonlinear response.

---

# 7. Rank‑7 transition and hysteresis in `g`

High-precision coexistence:

\[
\boxed{
g_c^{(7)}
=
3.718344898120381\ldots
}
\]

At coexistence localized profile:

\[
p_{\max}\approx0.836365.
\]

Full tangent Hessian at localized coexistence:
- positive,
- minimum ~`2.599064`.

Uniform Hessian at same `g`:
- also positive,
- minimum ~`3.290959`.

Hence first-order coexistence / metastability.

Refined landmarks:

\[
g_{\rm fold}\approx3.515644716839593,
\]

\[
g_c\approx3.718344898120381,
\]

\[
g_{\rm spin}\approx5.123427551398618.
\]

At fold:

\[
p_{\max}\approx0.6548232586.
\]

Metastable widths:

\[
g_c-g_{\rm fold}\approx0.2027001813,
\]

\[
g_{\rm spin}-g_c\approx1.4050826533,
\]

total hysteresis window:

\[
\approx1.6077828346.
\]

**Status:** rank7 branch/fold/coexistence is conditional/numerical except specific local repo certificates. Nie jest physical phase theorem.

---

# 8. Dual angular geometry — nested symmetry breaking

Fixed dual radius `r = ||theta||`, maximize:

\[
K(\theta)=\log\left[\frac1{12}\sum_i e^{\theta\cdot X_i}\right].
\]

Small `r`: angular maximizer = pure top `k=6` direction; gives 2 branches (even/odd).

At larger `r`: transition to 12 reflection-stabilized branches.

Group/orbit hierarchy:
- uniform stabilizer order 24 -> orbit size 1,
- pure `k6` stabilizer order 12 -> orbit size 2,
- full localized angular branch stabilizer order 2 -> orbit size 12.

Symbolically:

\[
1\to2\to12.
\]

Do not interpret literally as “bits” without additional model semantics.

## 8.1 Corrected pure-k6 angular instability

Earlier naive guess “k5 destabilizes k6 first” was corrected.

The first angular instability is controlled by `k3` due resonance:

\[
2k_3=k_6 \pmod{12}.
\]

Exact equation:

\[
\lambda_3(1+\tanh x)
=
\lambda_6\frac{\tanh x}{x},
\qquad
x=\sqrt{\lambda_6/12}\,r.
\]

Solution:

\[
r_{\rm spin}^{2\text{-branch}}
\approx0.4142113229119465.
\]

## 8.2 Angular first-order hysteresis

12-branch landmarks:

\[
r_{\rm fold}^{12}\approx0.3463027188406826,
\]

\[
r_{\rm coex}^{2\leftrightarrow12}
\approx0.36455557012840994,
\]

\[
r_{\rm spin}^{2}
\approx0.4142113229119465.
\]

Angular hysteresis width:

\[
\approx0.0679086041.
\]

At angular coexistence:
- 12-branch `pmax ~0.12174878`,
- pure-k6 `pmax ~0.09664001`,
- angular saddle `pmax ~0.11179786`,
- log-mgf barrier ~`1.29029e-5`.

Numerical fold transversality:
- `w^T F_r ~ -0.0059852`,
- `w^T F_zz[v,v] ~ 0.0317658`,
- coexistence crossing slope ~`0.001466`.

---

# 9. Cubic/quartic resonance structure and exact 12 phase locks

Po reflection reduction:

\[
z_3=r_3e^{i\phi_3},
\quad
z_4=r_4e^{i\phi_4},
\quad
z_5=r_5e^{i\phi_5},
\quad
z_6=s\,r_6,\quad s=\pm1.
\]

Cubic phase-sensitive resonances:

\[
3+3+6=12
\Rightarrow
s\,r_3^2r_6\cos(2\phi_3),
\]

\[
3+4+5=12
\Rightarrow
r_3r_4r_5\cos(\phi_3+\phi_4+\phi_5),
\]

\[
4+4+4=12
\Rightarrow
r_4^3\cos(3\phi_4).
\]

Cubic phase locks:

For `s=+`:
\[
\phi_3\in\{0,\pi\},
\quad
\phi_4=\frac{2\pi m}{3},
\quad
\phi_5=-\phi_3-\phi_4,
\quad m=0,1,2.
\]

6 solutions.

For `s=-`:
\[
\phi_3\in\{\pi/2,3\pi/2\},
\]
same `phi4`, and same sum rule.

6 more.

Total:

\[
\boxed{12\ \text{phase-locked states}}
\]

forming one D12 orbit; reflection stabilizer order 2.

Quartic phase resonances include:
- `cos 4phi3`,
- `cos(phi3-2phi4+phi5)`,
- `s cos(-phi3+phi4+phi5)`,
- `cos(phi3-3phi5)`,
- `s cos(phi4-2phi5)`.

At the cubic locked phases these quartic resonances are simultaneously optimized. Thus cubic+quartic fixed-positive-amplitude phase landscape has the same 12 locks.

---

# 10. Fourth cumulant almost reproduces full angular transition

Quartic truncation vs full log-mgf:

| landmark | quartic | full | relative error |
|---|---:|---:|---:|
| fold | 0.34660000115 | 0.34630271884 | +0.0858% |
| coexistence | 0.36428264569 | 0.36455557013 | -0.0749% |
| spinodal | 0.40987802444 | 0.41421132291 | -1.046% |

Interpretation:
- `kappa2`: spectral dominance,
- `kappa3`: creates resonant 12 phase orientations,
- `kappa4`: stiffens lock and generates first-order/metastable angular landscape,
- order 5+ corrections are small near angular transition.

At full angular coexistence representative amplitudes:

\[
|z_3|\approx0.1131879146,
\quad
|z_4|\approx0.1698528641,
\quad
|z_5|\approx0.2269339093,
\]

\[
z_6\approx-0.3380663037.
\]

Representative locked phases:

\[
\phi_3=-\pi/2,
\quad
\phi_4=2\pi/3,
\quad
\phi_5=-\pi/6.
\]

---

# 11. Morse census on phase 3-torus — major pre-Discord result

At fixed coexistence amplitudes and one sign of `z6`:

For quartic `K4`:
- exactly `60` distinct stationary roots found in 4000 random root solves,
- `6` maxima,
- `42` saddles,
- `12` minima,
- min abs nonzero Hessian eigenvalue ~`5.97884348e-5`.

For full:

- same census: `60 = 6+42+12`,
- min abs Hessian eigenvalue ~`5.29230377e-5`.

Matching quartic roots to full:
- all 60 matched,
- zero Morse-index changes,
- median phase displacement ~`2.96e-14` rad,
- max displacement ~`0.031890246` rad (~1.83 deg).

Both signs of `z6` -> total 12 angular maxima.

Sobol 65,536 phase points:

\[
\max|K_{\rm full}-K_4|
\approx2.1362\times10^{-6},
\]

\[
\max\|\nabla K_{\rm full}-\nabla K_4\|
\approx3.5722\times10^{-6},
\]

\[
\max\|H_{\rm full}-H_4\|_{\rm op}
\approx7.3461\times10^{-6}.
\]

Quartic stationary Hessian margin / correction:

\[
\frac{5.9788e-5}{7.3461e-6}
\approx8.14.
\]

Outside radius `0.08` critical neighborhoods:

sampled
\[
\min\|\nabla K_4\|
\approx1.1390e-5,
\]

gradient margin / correction ~`3.19`.

**Prospective theorem route:** interval-certify remainder `R_{>=5}`, complement gradient gap, and Krawczyk boxes around 60 roots. Then full phase topology would equal quartic topology at these amplitudes.

---

# 12. Trivial-stabilizer competitor audit

Repo ST344 proves symmetry alone cannot exclude generic orbit-24 stationary/minimizer candidates.

Our numerical scan across 18 radii `0.37 ... 13` found:
- 95 stationary orbit-24 competitors,
- `0` local angular maxima among found competitors.

Best gap to global angular max ~`0.0004815` at `r~0.456`.
Gap grows strongly:
- ~`0.518` near `r~3.6`,
- >`1` around `r~5`,
- >`4.5` near `r=13`.

**Status:** strong numerical evidence, not exhaustion theorem.

---

# 13. Radial nucleation vs angular branch selection

At `g=g_c^(7)`:

- angular branch selection occurs early around `r ~ 0.365`,
- main radial nucleation barrier around:
  \[
  r_\theta\approx1.81248,
  \]
- barrier:
  \[
  \Phi_{\rm barrier}\approx0.04655546,
  \]
- localized coexistence minimum:
  \[
  ||\theta||\approx3.537.
  \]

Important coordinate distinction:
- these are **theta-dual radii**;
- `h=X theta` has another norm; localized `||h||` was roughly `5.23`.
Do not mix them.

Interpretation inside conditional model: orientation/which branch can be decided long before strong concentration/localization.

---

# 14. Near-one-vertex feature ray

Full rank7 transition dual direction is almost aligned with one feature vector:

\[
\cos\angle(\theta,X_i)
\approx0.99992051,
\]

angle ~`0.72 deg`.

Restricting to ray:

\[
\theta=tX_0
\]

gives critical `g` only ~`0.0143%` above full 7D optimum.

This supports a strong radial core but **not** exact 1D reduction.

---

# 15. Repo guardrails relevant before Discord cutoff

Need preserve these exactly:

- ST344: open dense simplex stratum has trivial D12 stabilizer; symmetry does not force reflection-fixed minimizer.
- ST359: averaging by reflections can strictly raise `V4`; one-step rearrangement proof route fails.
- ST360: rank7 global uniqueness neither proved nor refuted.
- ST361: one positive rank7 reflection-even stationary root at `g=4` interval-certified; trivial-stabilizer competitors remain open.
- ST389/ST390: for **full strict A** at `g=4`, exactly 12 global minima in one D12 orbit.
- ST391/ST410: cap/positive-weight proof cannot be transferred to rank7 signed mediator; rank7 global minimality remains open.
- ST458: present sign-aware/absolute-value relaxations fail; this is not failure of the model.
- ST293/ST459: no strict active gain source found.

---

# PART II — `FIN_Discord_Robustness_and_Operational_Identifiability` AS CUTOFF

# 16. What the cutoff report established

The report used canonical:

\[
C=\frac I{12}+\frac W{20}.
\]

It constructs an explicit separable stationary state with marginals `C`, while proving canonical equilibrium with this marginal cannot be classical-quantum (zero one-sided discord class).

Important results from the report:
- separable state as mixture of 22 product density matrices,
- rank exactly 132,
- local heralded preparation with success probability `1/2`,
- entanglement not necessary,
- but zero-discord/CQ description excluded for canonical interaction,
- explicit constructed-state distance:
  \[
  D_{\rm tr}(R_{\rm sep},CQ_A)>1/15400,
  \]
- universal conservative lower bound for canonical equilibria:
  \[
  \sim3.45\times10^{-8}.
  \]

The report explicitly warned:
- preparation law/program is supplied,
- full high-order correlation law is not uniquely fixed by pairwise local geometry,
- no selector, units, apparatus, raw data, legacy role transfer, SM/GR or ToE closure.

This is the chronological cutoff requested by the user.

---

# PART III — BADANIA WYKONANE PO RAPORCIE DISCORD

# 17. New bridge: discord robustness is controlled by the same strict top spectral gap

From the report quantities:

\[
\Delta_L=0.1171091020573145
=
\frac{\lambda_6}{20},
\]

\[
\delta_L=0.00217878845336
=
\frac{\lambda_6-\lambda_5}{20}.
\]

Universal discord distance lower bound:

\[
d_0
=
\frac{\delta_L\Delta_L}{7392}
=
\frac{\lambda_6(\lambda_6-\lambda_5)}{2\,956\,800}
\approx3.45178516\times10^{-8}.
\]

Conditional active-selection linear growth:

\[
r_k=-1+\frac{g\lambda_k}{12}.
\]

Therefore:

\[
r_6-r_5
=
\frac g{12}(\lambda_6-\lambda_5)
=
\frac{5g}{3}\delta_L.
\]

Dimensionless top gap ratio:

\[
\frac{\lambda_6-\lambda_5}{\lambda_6}
\approx0.0186048.
\]

**Interpretation:** ten sam strict top eigenvalue/gap controls:
1. weak isolation of top active mode,
2. leading active growth-rate splitting in the conditional model,
3. robustness scale appearing in discord lower bound.

This is a structural bridge, not proof that discord *causes* localization or vice versa.

---

# 18. Exact CRT / ferromagnetic reduction of phase-locked rank7

After phase locking, active rank7 field can be represented on:

\[
\mathbb Z_4\times\mathbb Z_3.
\]

Let:

\[
\alpha=\frac\pi2(j\bmod4),
\qquad
\beta=-\frac{2\pi}{3}(j\bmod3).
\]

Then:

\[
\boxed{
h=
J_3\cos\alpha
+
J_4\cos\beta
+
J_5\cos(\alpha-\beta)
+
J_6\cos2\alpha
}
\]

with all `J_i >= 0` in the relevant locked sector.

Natural mapping:
- `k3`: `cos alpha`,
- `k4`: `cos beta`,
- `k5`: `cos(alpha-beta)`,
- `k6`: `cos 2alpha`.

This identifies the locked 4D exponential family as a cooperative/ferromagnetic finite model.

Consequences (Ginibre/Griffiths-type):
- relevant means nonnegative,
- mixed covariances nonnegative,
- Jacobian of fixed-point map is entrywise nonnegative,
- Perron-Frobenius leading direction can be chosen positive.

Fixed-point map:

\[
T(s)=g\nabla\log Z(s).
\]

Thus it is isotone in the positive orthant.

---

# 19. Cooperative fixed-point structure and collective saddle

Natural stationary states at rank7 coexistence:

Barrier saddle approximately:

\[
s_{\rm sad}
\approx
(0.94096,\ 1.00144,\ 0.96211,\ 0.68641),
\]

localized:

\[
s_{\rm loc}
\approx
(1.81990,\ 1.91399,\ 1.91457,\ 1.36720).
\]

Componentwise:

\[
0<s_{\rm sad}<s_{\rm loc}.
\]

Jacobian:

\[
DT=g\,\operatorname{Cov}(C).
\]

Eigenvalues:
- uniform: all < 1,
- localized: all < 1,
- saddle:
  \[
  0.3422,\ 0.4181,\ 0.5910,\ 1.2491.
  \]

Exactly one unstable direction numerically.

PF eigenvector at saddle roughly:

\[
v_{\rm PF}\approx(0.5286,0.5371,0.5348,0.3821).
\]

Alignment with actual displacement:

\[
\cos\angle(v_{\rm PF},s_{\rm loc}-s_{\rm sad})
\approx0.999554858.
\]

Interpretation: escape from barrier is collective positive activation of all active harmonics, not sign-frustrated single-mode motion.

---

# 20. Curvature strategy — target global MorseIndex <= 1

For natural 4D amplitude coordinates:

\[
H=
\frac1g I-\operatorname{Cov}_{p_s}(C).
\]

At coexistence `g=g_c` we want:

\[
\lambda_2(\operatorname{Cov}C)
<
\frac1{g_c}.
\]

Then Hessian has at most one negative eigenvalue globally:

\[
\boxed{\operatorname{MorseIndex}\le1}.
\]

Numerical critical level:

\[
1/g_c
\approx0.26893685965.
\]

A candidate global ceiling emerged:

\[
\boxed{
\sigma_*
=
\frac{
2\lambda_3(\lambda_4+\lambda_5)-\lambda_4\lambda_5
}{
24\lambda_3
}
}
\]

\[
\boxed{
\sigma_*
\approx0.26744324422884.
}
\]

Margin:

\[
1/g_c-\sigma_*
\approx0.00149361542
\]

(~`0.5554%`).

Associated:

\[
t_*^2=
\frac{(2\lambda_3-\lambda_4)(2\lambda_3-\lambda_5)}
{4\lambda_3^2},
\]

\[
t_*\approx0.4264779295.
\]

Natural `s3` coordinate on dangerous limit:

\[
s_3^*\approx0.796819388.
\]

---

# 21. Dangerous slice `s4=s5=0` — exact two-block structure

On:

\[
s_4=s_5=0,
\]

covariance splits into:
- `(k3,k6)` block,
- `(k4,k5)` block.

With:
\[
x=E[\cos\alpha],
\quad
y=E[\cos2\alpha],
\]

domain:

\[
0\le x\le1,
\qquad
x^2\le y\le1.
\]

Blocks:

\[
C_{36}=
\begin{pmatrix}
\frac{\lambda_3}{6}\left(\frac{1+y}{2}-x^2\right)
&
\sqrt{\frac{\lambda_3\lambda_6}{72}}x(1-y)
\\
\cdot&
\frac{\lambda_6}{12}(1-y^2)
\end{pmatrix},
\]

\[
C_{45}=
\begin{pmatrix}
\lambda_4/12&
\sqrt{\lambda_4\lambda_5}\,x/12\\
\cdot&\lambda_5/12
\end{pmatrix}.
\]

For `x<=t*`, `lambda_max(C45)<=sigma*`.

For `x>=t*`, algebraic PSD test:

\[
\sigma_*I-C_{36}\succeq0
\]

on whole `x^2<=y<=1`.

Smallest normalized determinant margin observed/algebraically reduced ~`0.02578`.

Equality only at:

\[
x=t_*,
\qquad
y=1,
\]

i.e.

\[
s_6\to\infty,
\quad
s_3=s_3^*,
\quad
s_4=s_5=0.
\]

This is a core exact-face result.

---

# 22. Boundary classification in compactified positive orthant

As fields -> infinity:

- `s3 -> inf`: support `{0,4,8}`, affine covariance rank 1 -> `lambda2=0`.
- `s4 -> inf`: support `{0,3,6,9}`, rank 2, ceiling `lambda6/12 ~0.19518`.
- `s5 -> inf`: single state -> `lambda2=0`.
- `s6 -> inf`: support `{0,2,4,6,8,10}`, rank 3 and can reach `sigma*`.
- multiple positive infinite fields collapse support further, usually rank <=1 / vertex.

Thus only `s6->inf` is dangerous for second curvature.

Broad 200k interior sample found no `lambda2 > sigma*`; max sampled ~`0.26150`.

Stationary searches for `grad lambda2=0` found no strictly interior candidate class; roots approached boundary faces.

**Status:** strong evidence, not yet global proof.

---

# 23. Important falsification: transverse monotonicity is false globally

Hypothesis:

> turning on `s4` or `s5` from `s4=s5=0` always lowers `lambda2`.

**[REFUTED].**

There are parts of the face where:

\[
\partial_{s_4^+}\lambda_2>0,
\quad
\partial_{s_5^+}\lambda_2>0.
\]

Max observed one-sided derivatives ~`0.08`.

However these regions are far below `sigma*`.

At balanced candidate `s3=s3*`, large `s6`, transverse derivatives are negative:

\[
d\lambda_2/ds_4\approx-0.03476,
\]

\[
d\lambda_2/ds_5\approx-0.03170,
\]

with negative quadratic coefficients too.

Thus candidate is a local ridge maximum on positive-orthant boundary, but no global simple monotonicity theorem.

---

# 24. Characteristic-polynomial route

For 4D covariance `M`, define:

\[
P(t)=\det(tI-M).
\]

At `t=sigma*`, target `lambda2<=sigma*` can be studied by sign changes:

\[
P(\sigma_*),P'(\sigma_*),P''(\sigma_*),P'''(\sigma_*),P''''.
\]

100k audit showed only sign patterns with <=1 relevant variation.

Potential proof chain:

\[
P\ge0\Rightarrow P'\ge0,
\]

\[
P'\ge0\Rightarrow P''\ge0,
\]

\[
P'''>0.
\]

## 24.1 Global `P''' > 0` via enclosing ball

Since:

\[
P'''(\sigma_*)=24\sigma_*-6\operatorname{tr}M,
\]

bound trace by minimum enclosing ball of 12 feature points.

Found center approximately:

\[
c=(0,\ 0.270913487,\ 0,\ 0.116305338),
\]

radius squared:

\[
R_{\rm enc}^2\approx0.927873534.
\]

\[
4\sigma_*\approx1.069772977.
\]

Thus:

\[
P'''(\sigma_*)
\ge
6(4\sigma_*-R^2)
\approx0.851397>0.
\]

Closed-form radius was derived.

## 24.2 Conditional `P' >=0 => P''>0`

Numerical constrained minimization:

\[
\min P''\ |\ P'\ge0
\approx0.035864947>0.
\]

Strong but initially numerical.

Hardest remaining implication in this route was:

\[
P\ge0\Rightarrow P'\ge0.
\]

This motivated Ising/probability reduction.

---

# 25. Exact two-spin Ising reduction on dangerous boundary `s6 -> inf`

On six surviving labels introduce binary:

\[
A\in\{\pm1\},
\quad
Y\in\{\pm1\},
\quad
B=\frac14+\frac34Y.
\]

Observables:

\[
k3\sim A,
\]

\[
k4\sim Y,
\]

\[
k5\sim\frac14A+\frac34 AY.
\]

Distribution:

\[
p(A,Y)
\propto
e^{H_AA+H_YY+KAY}
\]

with:

\[
H_A=J_3+\frac{J_5}{4},
\]

\[
H_Y=-\frac12\log2+\frac{3J_4}{4},
\]

\[
K=\frac{3J_5}{4}.
\]

Physical ferromagnetic domain:

\[
K\ge0,
\]

\[
H_Y\ge-\frac12\log2,
\]

\[
H_A\ge K/3.
\]

Partition function:

\[
Z=
2[
e^K\cosh(H_A+H_Y)
+
e^{-K}\cosh(H_A-H_Y)
].
\]

**Critical fact:** if these ferromagnetic constraints are relaxed, `lambda2` can exceed `sigma*`.
Therefore curvature ceiling is not pure geometry of support; it uses sign/coupling structure.

---

# 26. Semialgebraic probability coordinates for Ising boundary

Let:

\[
p_1=p_{++},\quad
p_2=p_{+-},\quad
p_3=p_{-+},\quad
p_4=p_{--}.
\]

Physical constraints become polynomial:

\[
p_1p_4-p_2p_3\ge0,
\]

\[
4p_1p_3-p_2p_4\ge0,
\]

\[
p_1p_2^2-p_3p_4^2\ge0,
\]

\[
p_i\ge0,\qquad \sum p_i=1.
\]

Characteristic coefficients:

\[
P(t)=t^3-e_1t^2+e_2t-e_3.
\]

`e1` = weighted pair-distance quadratic.

`e2` = cubic sum over tetrahedral face areas.

\[
e_3=
\frac{3\lambda_3\lambda_4\lambda_5}{8}
p_1p_2p_3p_4.
\]

Double-root candidate probabilities:

\[
p_{++}=\frac{1+t_*}{6},
\]

\[
p_{+-}=\frac{1+t_*}{3},
\]

\[
p_{-+}=\frac{1-t_*}{6},
\]

\[
p_{--}=\frac{1-t_*}{3}.
\]

At this point:

\[
P(\sigma_*)=0,
\qquad
P'(\sigma_*)=0,
\qquad
P''(\sigma_*)>0.
\]

---

# 27. Bernstein cover of Ising semialgebraic domain

Adaptive Bernstein cover on probability simplex:

- all boxes either:
  - excluded by physical constraints,
  - excluded by `P<0`,
  - or certified `P'>=0`,
- unresolved boxes collapse only around exact double-root.

At width ~`1/32768`:
- ~10 unresolved boxes,
- entire hull within ~`9e-5` of exact double-root.

No remote problematic region.

A second remote `P=0` component was numerically separated:
- nearest other `P=0` point distance ~`0.90568` in field coordinates,
- there `P' ~0.041995`.
- numerical minimum `P'` on remote component ~`0.04032345`, with KKT:
  \[
  \nabla P'\approx13.5998\,\nabla P.
  \]

Local cone around double-root:
- physical directions have strictly negative `P` second variation,
- full Hessian has an inadmissible direction requiring signs outside physical cone.

A local box `|X-X*|<=0.002`, `Y-1,Z-1 in [0,0.002]` was sufficient numerically/algebraically for supplied coefficients.

**Status:** boundary Ising theorem is essentially computer-assisted for supplied decimals; formal outward interval lift of strict eigenvalues still desirable.

---

# 28. Third covariance eigenvalue has lower ceiling

Found independent candidate ceiling:

\[
\boxed{
\sigma_3^*\approx0.2002057878
}
\]

with:

\[
q_3^*=
\frac{
3\lambda_4-\lambda_5+
\sqrt{(\lambda_5-\lambda_4)(16\lambda_3-9\lambda_4+\lambda_5)}
}{
6\lambda_4
}.
\]

Numerically:

\[
q_3^*\approx0.414684410.
\]

Then:

\[
\sigma_3^*
=
\frac{3\lambda_4}{8}q_3^*(1-q_3^*).
\]

Gap:

\[
\sigma_*-\sigma_3^*
\approx0.06723746.
\]

Numerical full fixed-`s6` optimization returns to `s3=s5=0` and finite `s4`.

Interpretation: even if second curvature approaches `sigma*`, the third curvature appears substantially lower. This supports effective single unstable direction geometry.

---

# 29. Exact parity decomposition of full 4D covariance

Split by `k6=+/-`.

Let:

\[
q=P(k6=+).
\]

Then exact law of total covariance:

\[
\boxed{
M
=
qC_+
+
(1-q)C_-
+
q(1-q)\Delta\mu\Delta\mu^T.
}
\]

Define:

\[
W=qC_++(1-q)C_-,
\]

\[
bb^T=q(1-q)\Delta\mu\Delta\mu^T.
\]

Thus:

\[
M=W+bb^T.
\]

Conditional structure:
- `C_+`: exact 3D Ising `(k3,k4,k5)` family.
- `C_-`: `k3` disappears; only 2D `(k4,k5)` sector.
- for fixed `(s3,s4,s5)`, changing `s6` changes only mixing weight `q`, not conditional covariance matrices.

This is central to all later proof reductions.

---

# 30. Exact first-order parity mixing away from double-root

At boundary double-root `q=1`, `C_+` has double eigenvalue `sigma*`.

Let:

\[
q=1-\varepsilon.
\]

Degenerate perturbation gives two first-order splits:

\[
\Delta\lambda_{k3}
=
\frac{\lambda_3(2t_*^2-1)}6
\varepsilon
\approx
-0.2079853448\,\varepsilon,
\]

\[
\Delta\lambda_{45}
=
-\frac{
\lambda_4\lambda_5t_*^2
}{
6\sqrt{(\lambda_5-\lambda_4)^2+4\lambda_4\lambda_5t_*^2}
}
\varepsilon
\approx
-0.0798064760\,\varepsilon.
\]

Thus both split **downward**.

After reoptimizing active `t` for each `q`:

\[
\boxed{
\lambda_{2,\max}(q)
=
\sigma_*
-
0.1312828584(1-q)
+
O((1-q)^2).
}
\]

So entering finite interior from the saturating boundary opens a linear gap.

---

# 31. Intraparity part `W` — Weyl closure architecture

Need prove:

\[
\lambda_2(W)\le\sigma_*.
\]

Exact parity weight bound at `s6=0`:

\[
\boxed{q_0\ge1/2}.
\]

Equality at `J3=J5=0` (arbitrary nonnegative `J4`).

For `C_-` exact maximal eigenvalue:

\[
\boxed{
\sup\lambda_1(C_-)
=
\frac{3\lambda_4+\lambda_5}{32}
=
0.27804102562746\ldots
}
\]

only slightly above `sigma*`.

Weyl:

\[
\lambda_2(W)
\le
q\lambda_2(C_+)
+
(1-q)\lambda_1(C_-).
\]

Two cases:

### Case A
If:
\[
\lambda_1(C_-)\le\sigma_*,
\]
then boundary Ising gives `lambda2(C+)<=sigma*`, hence immediate closure.

### Case B
If:
\[
\lambda_1(C_-)\ge\sigma_*,
\]
dangerous `C_-` forces strong localization in `C_+`.

---

# 32. Dangerous `C_-` -> dominant mass in `C_+`

Parameterize `C_-` probabilities:

\[
p_0=1-u,
\quad
p_+=\frac{u+d}{2},
\quad
p_-=\frac{u-d}{2}.
\]

Covariance:

\[
C_-=
\frac18
\begin{pmatrix}
3\lambda_4u(1-u)
&
-\sqrt{3\lambda_4\lambda_5}(1-u)d
\\
\cdot&
\lambda_5(u-d^2)
\end{pmatrix}.
\]

One can show `lambda_min(C-)<sigma*`, so dangerous condition equivalent to:

\[
\det(\sigma_*I-C_-)\le0.
\]

Determinant is linear in `d^2`, yielding exact minimal asymmetry `d_min(u)`.

Dangerous `u` interval:

\[
u_\pm=
\frac12
\left(
1\pm
\sqrt{
1-\frac{32\sigma_*}{3\lambda_4+\lambda_5}
}
\right),
\]

numerically:

\[
u_-\approx0.4023835428,
\quad
u_+\approx0.5976164572.
\]

Direct 1D minimum of dominant `C_+` mass:

\[
u_*\approx0.5280356754,
\]

\[
d_*\approx0.4522658586,
\]

\[
\boxed{
p_{\rm dom,min}
\approx0.719071046876.
}
\]

Required for geometry closure only:

\[
p_{\rm req}\approx0.7088170056.
\]

Margin:

\[
\approx0.0102540413.
\]

---

# 33. Conservative interval-style dominant-mass certificate

Instead of exact minimum, sufficient coarse bounds on entire dangerous interval:

\[
d/u\ge0.8515,
\]

\[
\frac{4(1-u)^2}{u^2-d^2}\ge11.45.
\]

These imply lower bounds on energy gaps and finally:

\[
\boxed{
p_{\rm dom}(C_+)\ge0.7112098557.
}
\]

Still > required `0.7088170056`.

Tetrahedral geometric envelope `e2_max(alpha)` is exactly decreasing for `alpha>1/2`:

\[
\frac{d e_{2,\max}}{d\alpha}
=
2Bx(1-2\alpha-2x)<0.
\]

At conservative `alpha`:

\[
e_2(C_+)\le0.06508301845,
\]

so:

\[
\lambda_2(C_+)
\le
\sqrt{e_2}
\approx0.2551137363.
\]

Then:

\[
\lambda_2(C_+)+\sup\lambda_1(C_-)
\le
0.5331547619
<
2\sigma_*
=
0.5348864885.
\]

Final margin:

\[
\boxed{
0.0017317265.
}
\]

**Interpretation:** if `C_-` becomes curvature-dangerous, shared fields force `C_+` to become concentrated enough to be harmless. The two parity sectors cannot be dangerous simultaneously.

This effectively closes `lambda2(W)<=sigma*` **conditional on boundary-Ising lemma** for supplied spectrum, with a robust margin.

---

# 34. Rank-one parity mixing — inertia / resolvent criterion

For:

\[
M=W+bb^T,
\]

if `W` has exactly one eigenvalue above `sigma*`, define:

\[
A=\sigma_*I-W.
\]

Rank-one inertia criterion:

\[
\boxed{
S
=
1-b^TA^{-1}b.
}
\]

A second supercritical direction appears only if:

\[
S<0.
\]

Large search:
- 200k points,
- full DE boxes up to `s_i<=12`,
- no `S<0`.

Global numerical minimizer:

\[
(s_3,s_4,s_5,s_6)
\approx
(1.703819,0,0,0),
\]

\[
\boxed{
S_{\min}
\approx0.057549460989.
}
\]

This is a large margin.

At minimizer:
- top `W` eigenvector lives in `(k4,k5)`,
- `b` lives in `(k3,k6)`,
- overlap ~zero.

Spectral identity:

\[
S=
1+
\frac{|b_1|^2}{\lambda_1-\sigma_*}
-
\sum_{i\ge2}
\frac{|b_i|^2}{\sigma_*-\lambda_i}.
\]

Thus overlap of `b` with already-supercritical top direction is protective (increases `S`). Worst geometry is orthogonality — exactly the observed face.

---

# 35. Exact 1D resolvent certificate on `s4=s5=s6=0`

Let:

\[
r=e^{-\sqrt{\lambda_3/6}s_3}\in(0,1).
\]

Then:

\[
S(r)=
1-
\frac{(\lambda_3/6)r(1-r)}
{\sigma_*(1+r)-(\lambda_3/6)r^2}
-
\frac{\lambda_6}{3\sigma_*}\frac{r}{(1+r)^2}.
\]

Common denominator:

\[
S(r)=
\frac{N(r)}
{(1+r)^2[\sigma_*(1+r)-(\lambda_3/6)r^2]}.
\]

Numerator cubic for supplied spectrum:

\[
N(r)
=
0.2674432442
-
0.3052987580r
-
0.6321999017r^2
+
0.8948405102r^3.
\]

Split `[0,1]` into `[0,1/2]`, `[1/2,1]`.

Bernstein coefficients all positive:
- first interval min coefficient ~`0.068599`,
- second ~`0.024204`.

Denominator factor:

\[
D_0(r)=
\sigma_*(1+r)-(\lambda_3/6)r^2
\]

is concave and positive at both endpoints:
- `D0(0)=sigma*>0`,
- `D0(1)=2sigma* - lambda3/6 ~0.20798534>0`.

Therefore:

\[
\boxed{
S(r)>0\quad \forall r\in(0,1)
}
\]

on entire face.

Derivative numerator quintic has numerically/Sturm exactly one root in `(0,1)`, giving unique face minimum:

\[
s_3\approx1.70381898,
\quad
S_{\min}\approx0.057549461.
\]

---

# 36. Transverse gap away from resolvent face

Full minimization forced away from face:

| required `max(s4,s5,s6)` | min S | gap above global min |
|---:|---:|---:|
| 0.01 | 0.0576410 | 9.16e-5 |
| 0.05 | 0.0590913 | 0.001542 |
| 0.10 | 0.0618139 | 0.004264 |
| 0.25 | 0.0718217 | 0.014272 |
| 0.50 | 0.0920439 | 0.034494 |
| 1.00 | 0.1655809 | 0.108031 |

At face minimizer:
- `dS/ds6 ~ +0.129578`,
- quadratic coefficients:
  - `H44 ~2.20888`,
  - `H45 ~2.22280`,
  - `H55 ~1.99223`.

Full 2x2 block is indefinite, but negative direction requires opposite signs `ds4 ds5 < 0`, outside physical positive cone. On `ds4,ds5>=0`, quadratic form positive.

---

# 37. Schur elimination of `k6`: resolvent -> ordinary 3x3 curvature

Because `W` has zero `k6` row/column, exact Schur complement:

\[
S>0
\iff
\lambda_2(\widetilde M)<\sigma_*,
\]

where:

\[
\boxed{
\widetilde M
=
W_{345}
+
\frac{
b_{345}b_{345}^T
}{
1-b_6^2/\sigma_*
}.
}
\]

This removes resolvent singularities and converts last mechanism into ordinary 3x3 second-eigenvalue problem.

---

# 38. Entire extreme face `s4=s5=0`, arbitrary `s3,s6` — algebraic certificate

Define:

\[
t=\tanh(\sqrt{\lambda_3/6}\,s_3),
\]

\[
q=P(k6=+).
\]

Physical `s6>=0` constraint becomes, with:

\[
x=qt,
\]

exactly:

\[
\boxed{
x^2\le2q-1,
\qquad q\in[1/2,1].
}
\]

On this face `tilde M` splits:
- one `k3` channel,
- `(k4,k5)` block.

The `(k4,k5)` top eigenvalue depends only on `x`:

\[
B(x)=
\frac{
\lambda_4+\lambda_5+
\sqrt{
(\lambda_5-\lambda_4)^2+
4\lambda_4\lambda_5x^2
}
}{24}.
\]

By definition of `t*`:

\[
x\le t_*
\Rightarrow
B(x)\le\sigma_*.
\]

For `x>=t*` only `k3` channel must be controlled.

Critical q:

\[
q_{\rm crit}
=
1-\frac1c
\approx0.6574434785,
\]

where:

\[
c=\frac{\lambda_6}{3\sigma_*}.
\]

Smallest q allowing `x>=t*`:

\[
q_{\min}
=
\frac{1+t_*^2}{2}
\approx0.5909417122.
\]

Two intervals:

### Interval 1: `q_min <= q <= q_crit`
Worst `x^2=2q-1`.
For polynomial:
\[
(\sigma_*-g(q))[1-cq(1-q)]
\]
Bernstein coefficients:

\[
(0.010316,\ 0.010091,\ 0.012743,\ 0.017993).
\]

Strictly positive.

### Interval 2: `q_crit <= q <= 1`
Worst `x^2=t_*^2`.
For:
\[
(\sigma_*-f(q))[1-cq(1-q)]
\]
Bernstein coefficients:

\[
(0.017993,\ 0.030538,\ 0.037327,\ 0).
\]

Nonnegative; equality only at `q=1`.

Therefore for supplied decimal spectrum:

\[
\boxed{
\lambda_2(\widetilde M)
\le\sigma_*
\quad
\forall s_3\ge0,\ s_6\ge0,\ s_4=s_5=0.
}
\]

Equality only in compactified boundary:

\[
q=1,\quad x=t_*,
\]

i.e. `s6->inf` at known double-root.

**This is one of the strongest post-Discord results.**

---

# 39. Current remaining global gap for MorseIndex <= 1

After all reductions, the difficult part is no longer the entire 11D simplex.

We now have:

1. **Boundary/phase-locked Ising ceiling** essentially certified for supplied spectrum.
2. **Intraparity W** effectively closed via two-case dominant-mass argument, conditional on boundary lemma.
3. **Rank-one mixing extreme face** algebraically certified for all `s3,s6`.
4. Strong numerical evidence that full 4D minimum/maximum returns to `s4=s5=0`.
5. Local physical-cone curvature says moving into positive `s4,s5` increases safety near identified extrema.

Remaining theorem gap:

\[
\boxed{
\text{exclude an off-face point }s_4>0\text{ or }s_5>0
\text{ with }\lambda_2(\widetilde M)>\sigma_*.
}
\]

Equivalent acceptable route:
- interval/Bernstein cover of compactified off-face domain,
- with exact face theorem as boundary condition,
- plus small local tube handled by copositive Hessian/one-sided derivatives.

If completed:

\[
\lambda_2(\operatorname{Cov}C)\le\sigma_*<1/g_c
\]

globally, hence:

\[
\boxed{
\operatorname{MorseIndex}(H)\le1
}
\]

for entire conditional rank7 active-gain landscape at coexistence.

Important:
- `MorseIndex<=1` **does not by itself prove global minimality / unique orbit**.
- It is a landscape-topology theorem, not selector theorem.
- It does not source gain `g`.

---

# 40. Important failed routes / do-not-repeat list

## [REFUTED] Symmetry-only reduction
Repo ST344 already proves generic trivial-stabilizer strata. Reflection-fixed search is not exhaustive.

## [REFUTED] One-step reflection averaging
ST359 gives state where every reflection averaging raises `V4`.

## [REFUTED] “s4/s5 always lower lambda2”
False away from balanced ridge.

## [REFUTED] Global monotonicity of optimized envelope in J4/J5
There are secondary branches/dips/rises, though all remain far below `sigma*`.

## [REFUTED] `lambda1+lambda2 <= 2 sigma*`
Counterexample:
\[
\lambda_1+\lambda_2\approx0.71612>2\sigma_*\approx0.53489.
\]

## [REFUTED] Simple block eigenvalue-repulsion shortcut
Naive bound using only max top eigenvalues of `(36)` and `(45)` blocks is not globally sufficient after mixed blocks turn on.

## [REFUTED / TOO LOOSE] Fixed Courant-Fischer hyperplanes
Candidates such as `u=mean`, `u=field`, `u=ones`, `u=C0-mu` do not globally certify required restricted variance.

## [REFUTED / TOO LOOSE] Fixed 3x3 principal minor interlacing
No single deleted coordinate gives a universal `<=sigma*` ceiling.

## [REFUTED / TOO LOOSE] Simple Weyl
\[
\lambda_2(W)+||bb^T||
\le\sigma_*
\]
is false/too weak.

## [REFUTED] Global `Q>=0` shortcut
Candidate:
\[
P' - (3/\sigma_*)P \ge0
\]
is false globally, although violations occur where `P<0`.

## [REFUTED] Lambda2 maximum over q always at endpoints
For some low-curvature configurations `lambda2(M(q))` has interior maximum. These maxima are far below `sigma*`.

## [REFUTED] Arbitrary relaxation of shared conditional fields
If the two parity classes are allowed independent `s4/s5`, resolvent `S<0` counterexamples appear. Physical shared exponential-family coupling is essential.

## [REFUTED] Arbitrary q below physical q0
If q is allowed below:
\[
q_0=Z_+/(Z_++Z_-),
\]
counterexamples to resolvent positivity appear. Physical q constraint matters.

---

# 41. Strong conceptual synthesis

The current conditional rank7 landscape exhibits several nested mechanisms:

## 41.1 Representation / active dimension
Strict rank7 field uses 7D spectral subspace.

## 41.2 Nonlinear entropy completion
Full 11D probability state is nonlinear softmax image; inactive low Fourier modes are generated automatically.

## 41.3 Angular orientation
At small radius: 2-branch pure `k6`.

At intermediate radius: resonant `k3,k4,k5,k6` phase lock creates 12 orientations.

## 41.4 Radial localization
Strong concentration occurs later after nucleation barrier.

## 41.5 Cooperative amplitude growth
CRT representation gives nonnegative mixed covariance; saddle escape is PF collective.

## 41.6 Curvature protection
Second unstable direction appears to be globally blocked by:
- competition between `(k3,k6)` and `(k4,k5)` channels,
- ferromagnetic coupling constraints,
- parity-sector complementarity,
- rank-one update geometry.

The dangerous point is not where one channel is maximized alone, but where two candidate second-curvature channels are exactly balanced:

\[
\lambda_{36}=\lambda_{45}=\sigma_*.
\]

This is a maximin/balanced-channel phenomenon.

---

# 42. Ontological/physical interpretation guardrail

What these results may support **inside the conditional model**:
- symmetry gives equivalent branch orbit,
- active positive feedback could amplify an initially tiny fluctuation,
- nonlinear entropy completion fills omitted coordinates,
- resonances organize 12 phase orientations,
- one effective unstable direction can mediate nucleation.

What they do **not** establish:
- source of `g`,
- why nature chooses rank7 mediator,
- physical clock/units,
- actual microscopic event law,
- apparatus/measurement,
- realized branch selector,
- QW-2191 discharge,
- Standard Model,
- GR,
- physical spacetime,
- `L_total`,
- ToE.

User ontology to preserve:
- nadsoliton itself is primordial information in solitonic state;
- no separate lower “information substrate” should be inserted;
- whole relational configuration is the nadsoliton, not a hub node.

---

# 43. Recommended exact next programme for second agent

## Priority 1 — outward-rounded upstream lift
Replace all frozen decimal `lambda3..lambda6` used in Bernstein coefficients by strict outward intervals or exact upstream enclosures.

Targets:
- `sigma*`,
- `t*`,
- Ising semialgebraic coefficients,
- dominant-mass thresholds,
- extreme-face Bernstein coefficients.

Goal: upgrade current “algebraic for supplied decimal spectrum” to theorem-grade computer-assisted proof.

## Priority 2 — off-face `s4,s5` global cover
Prove:
\[
\lambda_2(\widetilde M)\le\sigma_*
\]
for positive `s4,s5`.

Best architecture:
1. compactify amplitudes (`x_i = e^{-s_i}` or probability/simplex coordinates),
2. use exact face theorem at `s4=s5=0`,
3. local tube: one-sided/copolypositive Hessian,
4. complement: interval/Bernstein boxes,
5. explicitly include physical coupling constraints — do not relax classes independently.

## Priority 3 — exact phase Morse perturbation theorem
Independently complete the fixed-amplitude phase theorem:
- exact/interval classify 60 quartic roots,
- certify full order>=5 C2 remainder bounds,
- prove same Morse census for full log-mgf.

This attacks the trivial-stabilizer phase issue from functional structure rather than symmetry.

## Priority 4 — after MorseIndex theorem
Combine:
- stationary census,
- index data,
- parity/phase topology,
- compactness/boundary exclusion
to build a Morse/Conley connection graph.

Do **not** claim global rank7 minimality merely from MorseIndex<=1.

## Priority 5 — source problem remains logically separate
Do not invert landscape regularity into provenance of active gain.
ST293/ST459 guardrail remains.

---

# 44. Key numerical constants — quick reference

Strict:
- `lambda3 = 1.9614068619764449`
- `lambda4 = 2.199568849333209`
- `lambda5 = 2.2986062720790956`
- `lambda6 = 2.3421820411463004`

Rank7 transition:
- `g_fold = 3.515644716839593`
- `g_c = 3.718344898120381`
- `g_spin = 5.123427551398618`
- `pmax(coexist) ~ 0.836365`

Angular:
- `r_fold = 0.3463027188406826`
- `r_coex = 0.36455557012840994`
- `r_spin(k6) = 0.4142113229119465`

Curvature:
- `1/g_c ~ 0.26893685965`
- `sigma* = 0.26744324422884`
- margin ~`0.00149361542`
- `t* = 0.4264779295`
- dangerous natural `s3* ~ 0.796819388`

Third curvature:
- `sigma3* ~ 0.2002057878`

Parity:
- `sup lambda1(C-) = 0.27804102562746`
- direct dangerous dominant mass `p0_min ~ 0.719071046876`
- conservative p0 bound `>=0.7112098557`
- required `p0 >=0.7088170056`

Resolvent:
- `S_min ~0.057549460989`
- at `s3 ~1.70381898`, `s4=s5=s6=0`

Parity first-order split:
- `-0.2079853448 epsilon`
- `-0.0798064760 epsilon`

Optimized interior gap:
- `sigma* - lambda2_max(q) ~ 0.1312828584(1-q)`

Discord bridge:
- `Delta_L = lambda6/20 ~0.1171091020573`
- `delta_L = (lambda6-lambda5)/20 ~0.00217878845336`
- universal discord floor ~`3.451785e-8`.

---

# 45. Key generated artifacts from this chat

This handoff is self-contained, but the following generated files contain the computational ledgers/results and may be useful if available to the receiving agent.

### Pre-Discord / structural
- `FIN___kanoniczna_kowariancja_pr_d_w_cyklicznych.csv`
- `FIN___dok_adny_rozk_ad_fluktuacji_kraw_dziowych_11_55.csv`
- `FIN___pasywny_model_walker_w__szum_i_t_umienie_dok_adnie_si__r_wnowa__.csv`
- `FIN___strict_even_odd_memory_summary.csv`
- `FIN___no-go_dla_aktywnego_gainu_z_pasywnej_strict_pami_ci.csv`
- `FIN___exact_spectral-budget_threshold_at_rank_7.csv`
- `FIN___numeryczny_pierwszy_globalny_competitor_dla_ka_dego_aktywnego_ranku.csv`
- `FIN___max-entropy_completion_lowers_the_localization_threshold.csv`
- `FIN___nonlinear_Fourier_closure_of_cumulative_top-rank_sectors.csv`
- `FIN_rank-7___refined_hysteresis_landmarks.csv`
- `FIN_rank-7___numerical_C2_perturbation_margins_for_order_5_.csv`
- `FIN_rank-7___census_Morse_a_na_3-torusie_faz_przy_coexistence.csv`
- `FIN_rank-7___continuation_of_quartic_critical_points_into_full_log-mgf.csv`

### Post-Discord / curvature proof campaign
- `FIN_rank-7___exact_CRT_ferromagnetic_representation.csv`
- `FIN_rank-7___cooperative_fixed-point_Jacobians.csv`
- `FIN_rank-7___closed-form_asymptotic_second-curvature_candidate.csv`
- `FIN_rank-7___compactified_positive-orthant_boundary_census.csv`
- `FIN_rank-7___exact_two-spin_Ising_representation_on_the_dangerous_boundary.csv`
- `FIN_rank-7___semialgebraic_ferromagnetic_probability_domain.csv`
- `FIN_rank-7___Bernstein_cover_collapses_onto_the_double-root.csv`
- `FIN_rank-7___exact_first-order_parity-mixing_split.csv`
- `FIN_rank-7___conditional-parity_curvature_bounds.csv`
- `FIN_rank-7___direct_one-dimensional_dominant-mass_bound.csv`
- `FIN_rank-7___conservative_interval-style_closure_of_the_intraparity_W_bound.csv`
- `FIN_rank-7___one-dimensional_resolvent_minimum.csv`
- `FIN_rank-7___exact-face_Bernstein_positivity_certificate_for_the_resolvent.csv`
- `FIN_rank-7___corrected_Bernstein_certificate_for_the_full_extreme_face.csv`
- `FIN_rank-7___exact-face_certificate_landmarks.csv`

Key figures:
- `fin_rank7_full_vs_quartic_landmarks.png`
- `fin_rank7_phase_lock_curvature.png`
- `fin_rank7_quartic_structural_stability_margins.png`
- `fin_rank7_trivial_stabilizer_gap.png`
- `fin_rank7_monotone_fixed_point_iteration.png`
- `fin_rank7_lambda2_slice_to_infinity.png`
- `fin_rank7_bernstein_double_root_collapse.png`
- `fin_rank7_parity_mixing_split.png`
- `fin_rank7_intraparity_gap_vs_q.png`
- `fin_rank7_resolvent_margin_1d.png`
- `fin_rank7_resolvent_transverse_gap.png`
- `fin_rank7_direct_dominant_mass_1d.png`
- `fin_rank7_optimized_interior_gap.png`.

---

# 46. Provenance / cutoff sources used to build this handoff

Chronological cutoff report in this conversation:
- `Wklejony tekst(20260912-224119).txt`
  - title: FIN: separowalna równowaga i konieczny dysonans kwantowy
  - this is the material corresponding to the user's `FIN_Discord_Robustness_and_Operational_Identifiability` checkpoint.

Repo guardrail source available in conversation:
- `AGENTS(1).md`.

Other immediate pre-cutoff source:
- `Wklejony tekst(20260912-114934).txt`
  - microscopic/pure-equivalent interaction, rank-free strict-source obstruction and correlation-floor context.

This handoff intentionally distinguishes:
1. prior repo theorem,
2. new chat-derived exact algebra,
3. computer-assisted results,
4. strong numerical evidence,
5. falsified proof routes.

---

# 47. One-paragraph handoff for the receiving agent

The main new research line after the discord report is not a claim of new FIN physics but a near-closure of the **conditional rank‑7 active-gain landscape curvature problem**. The strict top Fourier rank‑7 mediator admits an exact 7D dual and, after phase locking, a 4D ferromagnetic CRT exponential family. Its first-order localization transition occurs numerically at `g_c=3.718344898120381`, with a nested angular 2→12 branch transition accurately controlled by cubic/quartic resonances. A 60-point phase Morse census is numerically stable under the full-vs-quartic remainder. For the amplitude landscape, a candidate global second-covariance ceiling
`σ*=0.26744324422884 < 1/g_c`
was derived in closed form. The dangerous boundary reduces exactly to a two-spin Ising model; semialgebraic/Bernstein analysis collapses all possible equality to one double-root. Full covariance splits by `k6` parity into `M=W+bb^T`; `W` is effectively controlled by a boundary-Ising theorem plus a dominant-mass complementarity argument, while the rank-one update has a large resolvent margin and an algebraically certified entire `s4=s5=0` face. The only substantial remaining curvature gap is an **off-face global exclusion for positive `s4` or `s5`**, ideally by outward-rounded interval/Bernstein cover. If completed, it would give global `MorseIndex<=1` for the supplied rank‑7 active-gain model, but it would still not source `g`, select a physical orbit member, prove rank‑7 global minimality, or close FIN physically.
