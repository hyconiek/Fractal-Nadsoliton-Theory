


To wina dwóch bardzo specyficznych błędów w parserze Markdown na GitHubie (MathJax):

1. Błąd **"Extra open brace..."** przy `\sum_{i<j}` wynika z tego, że znak `<` jest interpretowany przez Markdown jako... początek znacznika HTML (np. `<div`). To "psuje" klamry LaTeXa. Rozwiązaniem jest użycie komendy `\lt` (less than) zamiast znaku `<`.
2. Zepsute **3 punkty z listą** wynikają z tego, że GitHub bardzo źle radzi sobie z blokami wielolinijkowych równań umieszczonymi wewnątrz punktowanych list (gubi wcięcia). Rozwiązaniem jest zapisanie ich jako tzw. *inline math* (w jednej linii, otoczone pojedynczym znakiem `$`).

Wróciłem do standardowych bloków `$$`, naprawiłem usterkę z `<` oraz zamieniłem problematyczną listę na w pełni bezpieczny zapis w jednej linii. Teraz wyrenderuje się bezbłędnie.

Oto ostateczny kod:

***

```markdown
# A11 Nadsoliton Single‑Kernel Core Lagrangian + Emergence Map (Strict Candidate Packet)

Status: `A11_EXECUTED_NADSOLITON_SINGLE_KERNEL_CORE_LAGRANGIAN_AND_EMERGENCE_PACKET_NO_FALSE_PASS`  
As of: `2026-03-13`

## Goal

Export one self-contained **core** Lagrangian of the ToE in the repo’s own language:

```text
nadsoliton (one fundamental object) + one kernel as the internal coupling law,
and a strictly disciplined map of how “light → matter → observer” is meant to emerge from it,
without explaining the nadsoliton in terms of external theories.
```

This packet is a **definition/packaging** object:

1. it rewrites the already exported canonical $12\times\Psi + \Phi$ action template (`QW-2163/2165/2166`) into a compact
   “one kernel” form on a typed $\mathbb{Z}_{12}$ carrier, and
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
   - full canonical $12\times\Psi + \Phi$ Lagrangian density template with explicit kernel-index mixing symbols $K_{i,j}$.
2. `QW-2165`
   - exhaustive canonical EoM confirming the structural content (locality in $x$, self-polynomials, Yukawa cross terms,
     bidirectional kernel mixing).
3. `QW-2166`
   - exhaustive canonical Hessian / linearized EoM (diagonal stencil and $\Psi$--$\Phi$ cross-couplings).
4. `QW-2190/QW-2191`
   - the strict $n=12$ ring scaffold and the strict-core uniqueness obstruction from degenerate 2D modes ($O(2)$ family).
5. `QW-2118/QW-2049`
   - the strict working kernel tuple and the $n=12$ distance-profile evaluation used in the kernel-mode lane.
6. `F329`
   - typed $\mathbb{Z}_{12}$ carrier + regular action on the 12-slot scaffold (to name the internal “ring” indices).

## Scope & hard limits (no false pass)

This packet:

1. exports the **nadsoliton core** only (no external-theory interpretation),
2. does **not** assert that all local coefficient families (`m2_psi*`, `g4_psi*`, `g6_psi*`, `gY*`) are already
   strict-derived from the kernel alone,
3. does **not** claim a strict-core $O(2)$ cut on `pair1` or a discharge of `QW-2191` (see `T166` and the diagonal
   frontier),
4. does **not** add “half-life / decay-time” terms into the Lagrangian (open-system decay is not a conservative
   Lagrangian ingredient on the current strict scope).

## What follows from research (repo-state, no false PASS)

In the simplest possible terms, **at the current repo state**:

1. The nadsoliton core is already *formally* written as a local field system:
   - $12\times\Psi + \Phi$ with kinetics, local polynomial potential, and $\Phi^2 \Psi_i^2$ coupling,
   - and with **index mixing** $K_{i,j} \Psi_i \Psi_j$ in the same Lagrangian density (`QW-2163`),
   - and its EoM and Hessian/linearization are exported as strict objects (`QW-2165`, `QW-2166`).
2. In the kernel-mode lane, the repo possesses one selected strict working kernel $K_{\mathrm{sg}}(d)$ along with its profile on the
   $n=12$ ring (`QW-2049`, `QW-2118`). This is the current, hardest **candidate** for "one internal coupling law".
3. "Light" (first emergence) in the internal sense of the theory is:
   - **linearization modes** (eigenmodes of the operator/Hessian) around the nadsoliton vacuum,
   - which is an object strictly defined by `QW-2166` (stencil structure + mixing + cross-sections).
4. In the same lane, there is a hard uniqueness obstruction:
   - degenerate 2D pairs generate a continuous family of $O(2)$ basis choices, and the kernel itself does not pick an axis canonically
     (`QW-2191`).
5. If one wants a strict "choice accelerator" (cutting $O(2)$ on `pair1`), then mathematically:
   - a translationally invariant host is isotropic on `pair1` and does not cut $O(2)$ (`N465`),
   - the diagonal/local sector cuts $O(2)$ if and only if it has a non-zero mode-2 defect $F_2(d)$ (`N466`),
   - but **on current exports** $F_2(d)$ for the canonical `D_local_residual` remains undefined (`N472/P431`),
     so one must not declare a strict-core $O(2)$ cut.
6. A new, hard reducer (`N474`) states additionally:
   - under vacuum stationarity and $v_{\psi k} \neq 0$, the Yukawa contribution vanishes from the diagonal Hessian entry,
     so Yukawa cannot "single-handedly" provide the missing $F_2(d)$ in the strict core.
7. "Half-life / decay-time" is not yet a strict-core Lagrangian object:
   - there is no explicit damping in a conservative Lagrangian,
   - $t_{1/2}$ would require an object like width/instability in an effective description (open system),
   - so on the current strict scope, this can only be a *downstream interpretation*, not a term in $\mathcal{L}_{\mathrm{core}}$.
   - pure mathematical note: if an exponential decay
     $X(t)=X_0 e^{-\lambda t}$ is assumed downstream, then by definition $X(t_{1/2})=X_0/2$ yields
     $t_{1/2}=\frac{\ln 2}{\lambda}$; this is not a new ToE constant, but just a recalculation
     from the "half" condition.

## 1) Typed internal carrier and distance ($\mathbb{Z}_{12}$)

Let:

$$
I_{12}:=\{0,1,\ldots,11\}, \qquad \mathbb{Z}_{12}:=(I_{12},+ \bmod 12).
$$

Define the **directed** $\mathbb{Z}_{12}$ distance/step:

$$
d(i,j):=(j-i)\bmod 12 \in \{0,1,\ldots,11\}.
$$

This matches the strict kernel-mode lane convention where distance classes $1 \dots 11$ are evaluated as a profile
(`QW-2118`).

## 2) One strict working kernel (internal coupling law)

Define the strict working kernel:

$$
K_{\mathrm{sg}}(d) = \frac{\cos(\omega d+\phi)}{1+\beta d^{\eta}}, \qquad (\omega,\phi,\beta,\eta)=(0.18575, 0.16250, 1.0, 1.8).
$$

This is the later-pipeline strict working kernel selected by the strict gate chain (`QW-2049`) and used in the
kernel-mode ring lane (`QW-2118`).

Kernel-split discipline reminder:
this packet does **not** claim that $K_{\mathrm{sg}}$ has already inherited every historical role of any retired legacy kernel;
it only uses $K_{\mathrm{sg}}$ as the **single internal coupling law** of the present strict nadsoliton core candidate.

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

$\Phi$ is not treated as an independent ontological substrate; it is an internal order/coherence projection used in the
canonical action template.

## 4) The ToE core Lagrangian (one-kernel form)

### 4.1 Internal mixing matrix generated by the one kernel

Define the kernel-index mixing coefficients:

$$
K_{ij}:= \begin{cases} 0, & i=j \\ K_{\mathrm{sg}}(d(i,j)), & i\neq j \end{cases}
$$

Then the **kernel mixing potential** is:

$$
V_{\mathrm{mix}}(\Psi) := \frac{1}{2} \sum_{i\neq j} K_{ij}\,\psi_i\,\psi_j.
$$

Equivalently (grouping symmetric pairs). Note the use of `\lt` to fix GitHub parsing:

$$
V_{\mathrm{mix}}(\Psi) = \frac{1}{2} \sum_{i \lt j} (K_{ij}+K_{ji})\,\psi_i\,\psi_j.
$$

This form is chosen so that the Euler–Lagrange equation for $\psi_k$ contains the symmetrized coefficient
$(K_{k,j}+K_{j,k})/2$ exactly as exported in `QW-2165`.

### 4.2 Local self-interaction potential (canonical L13 template)

Define:

$$
V_{\Phi}(\phi) := \frac{1}{2} m_{\phi}^2\,\phi^2 + \frac{1}{4} \lambda_{\phi}\,\phi^4,
$$

$$
V_{\Psi}(\Psi) := \sum_{i=0}^{11}\left( \frac{1}{2} m_{\psi i}^2\,\psi_i^2 + \frac{1}{4} g4_{\psi i}\,\psi_i^4 + \frac{1}{6} g6_{\psi i}\,\psi_i^6 \right),
$$

$$
V_{\mathrm{Y}}(\Psi,\phi) := \sum_{i=0}^{11} gY_i\,\phi^2\,\psi_i^2.
$$

These coefficient families appear in the strict canonical action/EoM/Hessian exports (`QW-2163/2165/2166`).
This packet does **not** claim they are already strict-derived from the kernel alone.

### 4.3 Core Lagrangian density

Write the nadsoliton core Lagrangian density as:

$$
\mathcal{L}_{\mathrm{core}} = \frac{1}{2}\,\partial_\mu \phi\,\partial^\mu \phi + \frac{1}{2}\sum_{i=0}^{11}\partial_\mu \psi_i\,\partial^\mu \psi_i - V_{\Phi}(\phi) - V_{\Psi}(\Psi) - V_{\mathrm{Y}}(\Psi,\phi) - V_{\mathrm{mix}}(\Psi).
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
\psi_i(x)\equiv v_{\psi i}, \qquad \phi(x)\equiv v_{\phi}.
$$

The stationarity conditions are the constant-vacuum specializations of the exported EoM (`QW-2165`).

### 5.2 Light (first emergence) = linearized propagating modes

Linearize:

$$
\psi_i = v_{\psi i}+\eta_i, \qquad \phi = v_\phi+\eta_\phi.
$$

The quadratic action defines a linear operator (Hessian) on the fluctuation vector

$$
\eta = (\eta_0,\ldots,\eta_{11},\eta_\phi),
$$

with the **strictly exported** stencil structure (`QW-2166`):

1. diagonal $\Psi$ entries carry: $3g4_{\psi i}v_{\psi i}^2 + 5g6_{\psi i}v_{\psi i}^4 + 2gY_i v_\phi^2 + m_{\psi i}^2$
2. off-diagonal $\Psi$--$\Psi$ entries carry the symmetrized kernel mixing coefficients: $\frac{K_{ij}+K_{ji}}{2}$
3. $\Psi$--$\Phi$ cross terms carry: $4gY_i \, v_\phi \, v_{\psi i}$

“Light” in the internal ToE sense is:

```text
the propagating linearized eigenmodes of this operator around the nadsoliton vacuum.
```

On the $\mathbb{Z}_{12}$ ring scaffold (`QW-2190`), some eigenmodes come in degenerate 2D pairs, producing the strict $O(2)$
basis-choice family obstruction (`QW-2191`).

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

1. `QW-2191` remains a real strict-core uniqueness obstruction on degenerate mode pairs.
2. The diagonal/local accelerator lane reduces strict `pair1` $O(2)$ cutting to one mode‑2 defect decision target
   (`T166`).
3. Do not claim strict-core selector closure unless a genuinely strict internal selector source is exported.

## 7) What A11 does not prove

`A11` does not prove:

1. ToE closure,
2. strict-derived fixing of all local coefficient families from the kernel,
3. strict-derived vacuum values $(v_{\psi i}, v_\phi)$,
4. strict-core discharge of `QW-2191`,
5. any “half-life” law as a conservative-Lagrangian ingredient.

It only exports the cleanest honest core: **one nadsoliton**, **one kernel**, one explicit Lagrangian, and a disciplined
emergence map.
```