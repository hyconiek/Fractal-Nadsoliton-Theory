# A11 Nadsoliton Single‑Kernel Core Lagrangian + Emergence Map (Strict Candidate Packet)

Status: `A11_EXECUTED_NADSOLITON_SINGLE_KERNEL_CORE_LAGRANGIAN_AND_EMERGENCE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

Export one self-contained **core** Lagrangian of the ToE in the repo’s own language:

```text
nadsoliton (one fundamental object) + one kernel as the internal coupling law,
and a strictly disciplined map of how “light → matter → observer” is meant to emerge from it,
without explaining the nadsoliton in terms of external theories.
```

This packet is a **definition/packaging** object:

1. it rewrites the already exported canonical `12×Psi + Phi` action template (`QW-2163/2165/2166`) into a compact
   “one kernel” form on a typed `Z_12` carrier, and
2. it states (as a program-level emergence map) which operations on this Lagrangian correspond to the repo’s
   intended emergence ladder.

It does **not** claim ToE closure, strict-core selector closure, or `QW-2191` discharge.

## Ontology discipline (AX9)

The nadsoliton is the primordial information of the universe in a solitonic state.
There is no independent informational substrate underneath it.

Preferred internal order remains:

```text
nadsoliton → light → matter → emergent observer.
```

## Strict-admissible sources reused

1. `QW-2163`
   - full canonical `12×Psi + Phi` Lagrangian density template with explicit kernel-index mixing symbols `K_{i,j}`.
2. `QW-2165`
   - exhaustive canonical EoM confirming the structural content (locality in `x`, self-polynomials, Yukawa cross terms,
     bidirectional kernel mixing).
3. `QW-2166`
   - exhaustive canonical Hessian / linearized EoM (diagonal stencil and `Psi–Phi` cross-couplings).
4. `QW-2190/QW-2191`
   - the strict `n=12` ring scaffold and the strict-core uniqueness obstruction from degenerate 2D modes (`O(2)` family).
5. `QW-2118/QW-2049`
   - the strict working kernel tuple and the `n=12` distance-profile evaluation used in the kernel-mode lane.
6. `F329`
   - typed `Z_12` carrier + regular action on the 12-slot scaffold (to name the internal “ring” indices).

## Scope & hard limits (no false pass)

This packet:

1. exports the **nadsoliton core** only (no external-theory interpretation),
2. does **not** assert that all local coefficient families (`m2_psi*`, `g4_psi*`, `g6_psi*`, `gY*`) are already
   strict-derived from the kernel alone,
3. does **not** claim any **global** discharge of `QW-2191` from kernel-alone data; however, the current repo does export
   lane-scoped `O(2) -> Z2` axis cuts (diagonal/local and Shannon element-order reference lanes) which canonicalize all
   degenerate Fourier pairs on `n=12` **up to residual sign** (no selector closure implied),
4. does **not** add “half-life / decay-time” terms into the Lagrangian (open-system decay is not a conservative
   Lagrangian ingredient on the current strict scope).

## Co wynika z badań (repo-state, bez fałszywego PASS)

W możliwie prostych słowach, **na obecnym stanie repo**:

1. Rdzeń nadsolitonu jest już *formalnie* zapisany jako lokalny układ pól:
   - `12×Psi + Phi` z kinetyką, lokalnym potencjałem wielomianowym i sprzężeniem `Phi^2 Psi_i^2`,
   - oraz z **mieszaniem indeksowym** `K_{i,j} Psi_i Psi_j` w tej samej gęstości Lagrangianu (`QW-2163`),
   - a jego EoM i Hessian/linearizacja są wyeksportowane jako obiekty strict (`QW-2165`, `QW-2166`).
2. W lane kernel-mode repo posiada jeden wybrany strict kernel roboczy `K_sg(d)` wraz z jego profilem na pierścieniu
   `n=12` (`QW-2049`, `QW-2118`). To jest aktualny, najtwardszy **kandydat** na “jedno prawo sprzężeń wewnętrznych”.
3. „Światło” (pierwsze wyłonienie) w sensie wewnętrznym teorii to:
   - **mody liniaryzacji** (eigenmody operatora/Hessianu) wokół próżni nadsolitonu,
   - co jest obiektem ściśle zdefiniowanym przez `QW-2166` (struktura stencila + mieszania + przekroje).
4. W tej samej lane istnieje twarda przeszkoda unikatowości:
   - zdegenerowane pary 2D generują ciągłą rodzinę wyboru bazy `O(2)` i sam kernel nie wybiera osi kanonicznie
     (`QW-2191`).
5. Jednak strict core ma już dwa jawne, niezależne **wewnętrzne** mechanizmy osiowej kanonizacji (axis-only) na `n=12`,
   które tną `O(2)` do residual `Z2` na wszystkich parach `pair_m (m=1..5)` (bez globalnego discharge i bez selector closure):
   - diagonal/local: strict-derived value instantiation profilu `D_local_residual` przez `F447 → P437` rozstrzyga niezerowość
     defektów `F_{2m}(d)` na wszystkich parach (`N485`) i pakuje scoped `O(2)->Z2` discharge (`N487`) oraz obiekt bazy
     `ModeIndexAssignment_canonical_local_diagonal_strict_derived_v1` (`F453/N492`),
   - Shannon element-order reference: defekt `F_{2m}(ord_{Z_12})` oraz obiektyw cross-entropy tną `O(2)` do residual `Z2`
     na wszystkich parach (`N480`, `N488`, `N496`, wykonane przez `F454` i opakowane przez `N500`).
6. Reduktor `N474` pozostaje spójny z powyższym: diagonalny wpis Hessianu ma formę Yukawa-free pod stacjonarnością,
   więc Yukawa nie jest “źródłem osi” w tym strict mechanizmie; w szczególności diagonal/local `O(2)`‑cut nie wymaga
   osobnej “Yukawa‑orientacji”.
7. (Update, `2026-03-15`) Repo eksportuje też pełny, liczbowy eigensystem *pełnej* macierzy Hessianu sektora `Psi`
   na diagonal/local lane:
   `H_psi := K_total + (m0^2 I + D_local_residual)` wraz z ortonormalną bazą własną i wartościami własnymi (`F459`).
   To jest strict-derived value instantiation (lane-scoped) i wspiera interpretację “light = linearized eigenmodes”
   bez promocji do host matching / ToE closure.
8. (Update, `2026-03-15`) Aby nie przemycać “fizycznej orientacji” przez arbitralny znak wektorów własnych,
   repo eksportuje sign‑gauge‑invariant rank‑one projektory spektralne `P_j := |v_j><v_j|` dla `H_psi` (`F460`) i pakuje
   gauge‑irrelewantność znaku jako theorem (`N504`).
9. (Update, `2026-03-15`) Na aktualnej instancji `H_psi` eigenmody są silnie zmieszane pomiędzy Fourier pair planes:
   sonda eksportuje profil wag `w_{j,label} := tr(Π_label P_j)` po labelach `{e0,pair1..pair5,e6}` (`P464`), a wniosek jest
   opakowany jako value‑instantiation theorem (`N505`). To blokuje jakikolwiek “cichy” skrót typu “pair-plane modes diagonalize H_psi”.
10. „Czas połowicznego rozpadu” nie jest jeszcze obiektem strict-core Lagrangianu:
   - w konserwatywnym Lagrangianie nie ma wprost tłumienia,
   - `t_{1/2}` wymagałby obiektu typu szerokość/niestabilność w opisie efektywnym (otwarty układ),
   - więc na obecnym strict scope to może być tylko *interpretacja downstream*, nie termin w `L_core`.
   - uwaga czysto matematyczna: jeśli downstream przyjmuje się wykładniczy zanik
     $X(t)=X_0 e^{-\lambda t}$, to z definicji $X(t_{1/2})=X_0/2$ daje
     $t_{1/2}=\frac{\ln 2}{\lambda}$; to nie jest nowa stała ToE, tylko przeliczenie
     z warunku „połowy”.

## 1) Typed internal carrier and distance (Z\_12)

Let:

$$
I_{12}:=\{0,1,\ldots,11\},
\qquad
\mathbb{Z}_{12}:=(I_{12},+ \bmod 12).
$$

Define the **directed** `Z_12` distance/step:

$$
d(i,j):=(j-i)\bmod 12 \in \{0,1,\ldots,11\}.
$$

This matches the strict kernel-mode lane convention where distance classes `1..11` are evaluated as a profile
(`QW-2118`).

## 2) One strict working kernel (internal coupling law)

Define the strict working kernel:

$$
K_{\mathrm{sg}}(d)
=
\frac{\cos(\omega d+\phi)}{1+\beta d^{\eta}},
\qquad
(\omega,\phi,\beta,\eta)=(0.18575,\ 0.16250,\ 1.0,\ 1.8).
$$

This is the later-pipeline strict working kernel selected by the strict gate chain (`QW-2049`) and used in the
kernel-mode ring lane (`QW-2118`).

Kernel-split discipline reminder:
this packet does **not** claim that `K_sg` has already inherited every historical role of any retired legacy kernel;
it only uses `K_sg` as the **single internal coupling law** of the present strict nadsoliton core candidate.

## 3) Nadsoliton core fields

Introduce:

1. a 12-component real field on spacetime (the nadsoliton carrier degrees of freedom):
   $$
   \Psi(x) = (\psi_0(x),\ldots,\psi_{11}(x))\in\mathbb{R}^{12},
   $$
2. one real scalar order / coherence field:
   $$
   \Phi(x)=\phi(x)\in\mathbb{R}.
   $$

`Phi` is not treated as an independent ontological substrate; it is an internal order/coherence projection used in the
canonical action template.

## 4) The ToE core Lagrangian (one-kernel form)

### 4.1 Internal mixing matrix generated by the one kernel

Define the kernel-index mixing coefficients:

$$
K_{ij}:=
\begin{cases}
0, & i=j,\\[2mm]
K_{\mathrm{sg}}(d(i,j)), & i\neq j.
\end{cases}
$$

Then the **kernel mixing potential** is:

$$
V_{\mathrm{mix}}(\Psi)
:=
\frac12 \sum_{i\neq j} K_{ij}\,\psi_i\,\psi_j.
$$

Equivalently (grouping symmetric pairs):

$$
V_{\mathrm{mix}}(\Psi)
=
\frac12 \sum_{i<j} (K_{ij}+K_{ji})\,\psi_i\,\psi_j.
$$

This form is chosen so that the Euler–Lagrange equation for $\psi_k$ contains the symmetrized coefficient
$(K_{k,j}+K_{j,k})/2$ exactly as exported in `QW-2165`.

### 4.2 Local self-interaction potential (canonical L13 template)

Define:

$$
V_{\Phi}(\phi)
  :=
  \frac12 m_{\phi}^2\,\phi^2 + \frac14 \lambda_{\phi}\,\phi^4,
$$

$$
V_{\Psi}(\Psi)
  :=
  \sum_{i=0}^{11}\left(
    \frac12 m_{\psi i}^2\,\psi_i^2
    +\frac14 g4_{\psi i}\,\psi_i^4
    +\frac16 g6_{\psi i}\,\psi_i^6
  \right),
$$

$$
V_{\mathrm{Y}}(\Psi,\phi)
  :=
  \sum_{i=0}^{11} gY_i\,\phi^2\,\psi_i^2.
$$

These coefficient families appear in the strict canonical action/EoM/Hessian exports (`QW-2163/2165/2166`).
This packet does **not** claim they are already strict-derived from the kernel alone.

### 4.3 Core Lagrangian density

Write the nadsoliton core Lagrangian density as:

$$
\mathcal{L}_{\mathrm{core}}
=
\frac12\,\partial_\mu \phi\,\partial^\mu \phi
+
\frac12\sum_{i=0}^{11}\partial_\mu \psi_i\,\partial^\mu \psi_i
-
V_{\Phi}(\phi)
-
V_{\Psi}(\Psi)
-
V_{\mathrm{Y}}(\Psi,\phi)
-
V_{\mathrm{mix}}(\Psi).
$$

On the strict L13 variational gate exports (`QW-2163/2165/2166`), one may read the same density on the 1D slice
with $\partial_\mu\partial^\mu \to d^2/dx^2$ (this is an execution lane, not a claim that the full Lorentzian package
is already closed).

## 5) What “emergence” means here (map from the Lagrangian)

This section is a disciplined **map** of operations on the core Lagrangian.
It is not a claim that every arrow is already theorem-level discharged.

### 5.1 Nadsoliton vacuum (background)

Choose a constant vacuum:

$$
\psi_i(x)\equiv v_{\psi i},
\qquad
\phi(x)\equiv v_{\phi}.
$$

The stationarity conditions are the constant-vacuum specializations of the exported EoM (`QW-2165`).

### 5.2 Light (first emergence) = linearized propagating modes

Linearize:

$$
\psi_i = v_{\psi i}+\eta_i,
\qquad
\phi = v_\phi+\eta_\phi.
$$

The quadratic action defines a linear operator (Hessian) on the fluctuation vector

$$
\eta = (\eta_0,\ldots,\eta_{11},\eta_\phi),
$$

with the **strictly exported** stencil structure (`QW-2166`):

1. diagonal `Psi` entries carry:
   $$
   3g4_{\psi i}v_{\psi i}^2 + 5g6_{\psi i}v_{\psi i}^4 + 2gY_i v_\phi^2 + m_{\psi i}^2,
   $$
2. off-diagonal `Psi–Psi` entries carry the symmetrized kernel mixing coefficients:
   $$
   \frac{K_{ij}+K_{ji}}{2},
   $$
3. `Psi–Phi` cross terms carry:
   $$
   4gY_i\,v_\phi\,v_{\psi i}.
   $$

“Light” in the internal ToE sense is:

```text
the propagating linearized eigenmodes of this operator around the nadsoliton vacuum.
```

On the `Z_12` ring scaffold (`QW-2190`), some eigenmodes come in degenerate 2D pairs, producing the strict `O(2)`
basis-choice family obstruction (`QW-2191`).

Update (`2026-03-15`): na obecnym repo state kernel-alone `QW-2191` pozostaje prawdziwe jako obstruction, ale w strict core
istnieją lane-scoped wewnętrzne składniki, które kanonizują osie w parach (do residual `Z2`) i umożliwiają jawne,
liczbowe pakiety eigensystemów (np. diagonal/local `F459`).

### 5.3 Matter (second emergence) = stable nonlinear excitations of the same fields

“Matter” is intended to arise as **nonlinear**, spatially localized and/or topologically stabilized excitations of the
same core fields, supported by:

1. the self-interaction nonlinearities in $V_{\Psi}$ and $V_{\mathrm{Y}}$,
2. the internal kernel-induced mixing structure in $V_{\mathrm{mix}}$,
3. the existence and stability analysis of nontrivial solutions to the Euler–Lagrange system.

This packet does not claim the full classification of such excitations is already strict-closed. It only fixes the
precise core object they must come from: the same $\mathcal{L}_{\mathrm{core}}$.

### 5.4 Emergent observer (third emergence) = self-referential closure patterns (not a new substrate)

An “observer” is not introduced as a new ontological layer.

It is intended to emerge as a late-stage, self-referential closure pattern supported by sufficiently complex and
stable matter-sector excitations of the same nadsoliton fields.

This is a program-level ontology statement; it is not promoted here into a theorem-level closure claim.

## 6) Current strict frontier reminders (what is still open)

1. `QW-2191` remains a real strict-core uniqueness obstruction on kernel-alone translation-invariant data.
2. The diagonal/local `T166` decision target is now discharged on a strict-derived value instantiation (`N482`), and the
   same lane exports an all-pairs `O(2)->Z2` cut on `n=12` (`N487`) plus explicit mode-index assignment basis (`F453`).
3. Do not claim strict-core selector closure / admissible `S_sel_int` unless a genuinely strict closure object is exported.

## 7) What A11 does not prove

`A11` does not prove:

1. ToE closure,
2. strict-derived fixing of all local coefficient families from the kernel,
3. strict-derived vacuum values $(v_{\psi i},v_\phi)$,
4. strict-core discharge of `QW-2191`,
5. any “half-life” law as a conservative-Lagrangian ingredient.

It only exports the cleanest honest core: **one nadsoliton**, **one kernel**, one explicit Lagrangian, and a disciplined
emergence map.
