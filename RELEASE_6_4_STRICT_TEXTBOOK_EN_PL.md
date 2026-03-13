# RELEASE 6.4 STRICT TEXTBOOK EDITION (EN + PL)

**Version:** 6.4.2  
**Date:** 2026-03-13  
**Branch:** `main`  
**Predecessor:** Release 6.3 — Global Selector Capture Edition

**Status discipline (no false pass):**
- This document is a **strict-only textbook projection** of the current repo state.
- It does **not** claim strict-core theta export, strict-core selector closure, global selector closure,
  `QW-2191` discharge, or ToE closure.
- Any use of extension-lane representatives (e.g. fixing selector slots) remains explicitly labeled
  `strict_extension_only` and must not be promoted into strict core.

---

## ENGLISH VERSION

## 0) Delta Since Release 6.3

Release 6.4 adds **new strict structural closures** around the `QW-2191` uniqueness obstruction:

1. the certified translation-invariant host operator is **provably isotropic** on `pair1 = span{c1,s1}` and cannot
   cut the `O(2)` family (`N465`, audited by `P424`),
2. a diagonal/local sector can cut `O(2)` on `pair1` **iff** its diagonal profile has a nonzero **mode‑2** Fourier
   coefficient `F2(d)` (`N466`, audited by `P425`),
3. for `n=12`, the relevant diagonal mode‑2 defect reduces to **six opposite-pair sums**
   `S_k := d_k + d_{k+6}` (exactly the six classes already exported by `R18`), hence the strict `pair1` diagonal
   `O(2)`-cut question reduces to one explicit checkable complex linear combination (`N467`, persisted as `P426`).
4. if `F2(d)≠0`, the induced `pair1` diagonal axis/eigenbasis is explicit and canonical:
   $\theta_*=\frac12\,\mathrm{atan2}(\mathrm{Im}\,F_2,\ \mathrm{Re}\,F_2)$ (`N468`, audited by `P427`),
5. the strict sigma-int FR-derived $\mathbb{Z}_2$ parity sign mask $b_k:=(-1)^k$ has **zero** mode‑2 defect
   $F_2(b)=0$, hence cannot cut `O(2)` on `pair1` by itself (`N469`, audited by `P428`),
6. any +3-shift-invariant diagonal profile ($d_{k+3}=d_k$ on `n=12`) has **zero** mode‑2 defect $F_2(d)=0$, hence
   cannot cut `O(2)` on `pair1` (`N470`, audited by `P429`).
7. for `n=12`, if a diagonal profile is `Aut(Z_12)`-invariant (constant on the 6 quotient orbits from `N455`), then
   its mode‑2 defect is forced **real** and reduces to one explicit orbit-linear combination; hence the diagonal `pair1`
   `O(2)`-cut question collapses to one real scalar check on six orbit values (`N471`, audited by `P430`).
8. on the current exported canonical FIN diagonal coefficient class (`R15`), `F2(d)` for the canonical local diagonal
   residual sector is **underdetermined**: both `F2(d)=0` and `F2(d)≠0` remain compatible with current exports (`N472`,
   audited by `P431`).
9. the repo currently exports **no numeric instantiation** of the six opposite‑pair sums (or equivalent diagonal values)
   sufficient to evaluate `F2(d)` directly (`P432`).
10. the mode‑2 defect can be recovered equivalently from the `R18` declared `pair1` entry system
   (`N473`, audited by `P433`), and an explicit evaluation harness is exported for when strict-derived opposite‑pair
   sums become available (`P434`).
11. under constant-vacuum stationarity (exported canonical EoM) and `vpsi_k ≠ 0`, the Yukawa combination cancels out of
   the canonical diagonal entry defining `d_k`, yielding a Yukawa-free diagonal expression (`N474`, toy-audited by
   `P435`).

This is a *reduction*, not a discharge: it makes the “physical accelerator of choice” claim checkable without
introducing new hidden slots, but it does not yet provide strict-derived diagonal coefficient values. In particular,
`N472/P431` show the canonical coefficient class does not decide `F2(d)` on current exports, and `P432` confirms no
value instantiation is currently exported to decide it by direct evaluation.

## 1) One-Page Strict Status (6.4)

### 1.1 What is now structurally sharp

1. `QW-2191` remains true: kernel alone leaves a continuous `O(2)` basis-choice family in each degenerate 2D mode
   subspace.
2. The host kernel route cannot solve this on `pair1`: the certified host operator
   `A_host := K_total + m0^2 I` is isotropic on `pair1` (`N465`).
3. Any strict `pair1` `O(2)` cut must come from **non-translation-invariant structure**, most naturally a diagonal/local
   sector with `F2(d) ≠ 0` (`N466`).
4. For the canonical FIN local diagonal residual sector on the 12-slot carrier, the full `pair1` diagonal `O(2)`-cut
   condition reduces to a 6-class mode‑2 defect expression (no hidden “choose a site” slot): `N467/P426`.
5. If such a diagonal cut exists (`F2(d)≠0`), the diagonalization axis on `pair1` is canonical and computed directly
   from `F2(d)` (`N468/P427`).
6. Two common “symmetry breaking” slogans are now closed as diagonal `pair1` `O(2)`-cut sources:
   - the FR parity sign mask $(-1)^k$ has `F2=0` (`N469/P428`),
   - +3 shift invariance ($d_{k+3}=d_k$) forces `F2=0` (`N470/P429`).
7. If the diagonal/local profile is additionally `Aut(Z_12)`-invariant (a quotient-safe, no-hidden-slot condition),
   then `Im(F2)=0` and the `pair1` diagonal `O(2)`-cut criterion reduces to one real scalar check on six orbit values
   (`N471/P430`).
8. On the current canonical FIN diagonal coefficient class (`R15`), `F2(d)` for the canonical local diagonal residual
   sector is **underdetermined** at the level of exports, hence no strict-core diagonal `pair1` `O(2)`-cut witness can
   be promoted yet (`N472/P431`).

### 1.2 What is still missing (the real frontier)

The repo still does **not** export:

1. any strict-derived diagonal/local coefficient/value instantiation deciding whether the canonical local diagonal
   sector has `F2(d) ≠ 0` (currently underdetermined on exports: `N472/P431`),
2. any strict-core canonical theta-supply ingredient upgrading candidate representatives into strict core (`T159`),
3. any global `QW-2191` discharge or ToE closure.

So Release 6.4 is a strict reduction step: it replaces vague “selector” rhetoric with an explicit `F2(d)` defect
question.

## 2) Strict Kernel (Operational Gate Kernel)

The current strict working kernel is:

$$
K_{\mathrm{sg}}(d)
=
\frac{\cos(\omega d+\phi)}{1+\beta d^{\eta}},
\qquad
(\omega,\phi,\beta,\eta)=(0.18575,\ 0.16250,\ 1.0,\ 1.8).
$$

Define:

$$
\Theta(d):=\omega d+\phi,
\qquad
D(d):=1+\beta d^{\eta} = 1+d^{1.8},
$$

so:

$$
K_{\mathrm{sg}}(d)=\frac{\cos(\Theta(d))}{D(d)}.
$$

## 3) `QW-2191` and the `pair1` `O(2)` family (context)

On the 12-octave ring (`QW-2190`), each degenerate 2D mode subspace admits a continuous basis rotation:

$$
(c_1,s_1)\ \mapsto\ (c_1',s_1')=(c_1,s_1)\,R(\theta),
\qquad
R(\theta)=\begin{pmatrix}\cos\theta & -\sin\theta\\ \sin\theta & \cos\theta\end{pmatrix}.
$$

`QW-2191` states that kernel-only data does not fix $\theta$ canonically; therefore any strict-core uniqueness must
export an additional internal selector source or an explicit symmetry-breaking premise.

Release 6.4 focuses on one precise “accelerator” subquestion:

```text
Can the certified host operator (or the canonical local diagonal residual sector) itself break O(2) on pair1?
```

## 4) New strict result A: host isotropy on `pair1` (N465 + P424)

Let the certified host operator be:

$$
A_{\mathrm{host}} := K_{\mathrm{total}} + m_0^2 I.
$$

### 4.1 Theorem (N465)

On the current strict host route, the restriction of $A_{\mathrm{host}}$ to `pair1 = span{c_1,s_1}` is scalar:

$$
\left.A_{\mathrm{host}}\right|_{\mathrm{pair1}} = \lambda_1\, I_2.
$$

Therefore it cannot supply an `O(2)` cut on `pair1`.

### 4.2 Numeric audit (P424)

`P424` computes the explicit `pair1` block in basis `(c1,s1)`:

$$
\left[A_{\mathrm{host}}\right]_{\{c_1,s_1\}}
=
\begin{pmatrix}
1.9197380969174072 & 5.144911710811878\times 10^{-17}\\
5.144911710811878\times 10^{-17} & 1.9197380969174067
\end{pmatrix},
$$

with anisotropy signature:

$$
\Delta_1(A_{\mathrm{host}}) = (a_1-d_1,\ b_1)
\approx
(4.44\times 10^{-16},\ 5.14\times 10^{-17}),
$$

consistent with isotropy at tolerance $10^{-12}$.

So any attempted strict `pair1` axis selection must come from additional non-translation-invariant structure.

## 5) New strict result B: exact diagonal `pair1` `O(2)`-cut criterion (N466 + P425)

Let $D=\mathrm{diag}(d_0,\ldots,d_{n-1})$ be a diagonal operator in the site basis on an $n$-ring.

Define the mode‑2 Fourier coefficient:

$$
F_2(d) := \sum_{i=0}^{n-1} d_i\, e^{i\frac{4\pi i}{n}} \in \mathbb{C}.
$$

### 5.1 Theorem (N466)

On `pair1`, the diagonal restriction is scalar **iff** $F_2(d)=0$.
Equivalently, the `pair1` anisotropy signature is:

$$
\Delta_1(D) = (a_1-d_1,\ b_1) =
\left(\frac{2}{n}\,\mathrm{Re}\,F_2(d),\ \frac{1}{n}\,\mathrm{Im}\,F_2(d)\right).
$$

So a diagonal/local sector breaks `O(2)` on `pair1` **iff** $F_2(d)\neq 0$.

### 5.2 Numeric audit (P425)

`P425` checks this criterion on explicit toy profiles for $n=12$:

- constant profile: $\Delta_1 \approx (3.33\times 10^{-16},\ 2.04\times 10^{-17})$ (no cut),
- mode‑1 profile: $\Delta_1 \approx (-2.40\times 10^{-16},\ -8.37\times 10^{-17})$ (no cut),
- mode‑2 cosine profile: $\Delta_1 \approx (1.0,\ 0)$ (cuts `O(2)`),
- mode‑2 sine profile: $\Delta_1 \approx (0,\ 0.5)$ (cuts `O(2)`).

## 6) New strict result C: canonical local diagonal mode‑2 defect reduces to 6 classes (N467 + P426)

Now specialize to the FIN canonical 12-slot carrier (indices `k=0..11`).

Let the diagonal profile of the **canonical local diagonal residual sector** be:

$$
d_k := (D_{\mathrm{loc}})_{kk},
\qquad k=0,\ldots,11.
$$

Define opposite-pair sums:

$$
S_k := d_k + d_{k+6},
\qquad k=0,\ldots,5.
$$

### 6.1 Theorem (N467)

For $n=12$:

$$
F_2(d)
\;=\;
\sum_{i=0}^{11} d_i\,e^{i\frac{4\pi i}{12}}
\;=\;
\sum_{i=0}^{11} d_i\,e^{i\frac{\pi i}{3}}
\;=\;
\sum_{k=0}^{5} (d_k+d_{k+6})\,e^{i\frac{\pi k}{3}}
\;=\;
\sum_{k=0}^{5} S_k\,e^{i\frac{\pi k}{3}}.
$$

Writing the six phase factors out yields the explicit reduction:

$$
\mathrm{Re}\,F_2
=
(S_0 - S_3) + \frac{1}{2}(S_1 - S_2 - S_4 + S_5),
$$

$$
\mathrm{Im}\,F_2
=
\frac{\sqrt{3}}{2}(S_1 + S_2 - S_4 - S_5).
$$

So the strict `pair1` diagonal `O(2)`-cut question for the canonical local diagonal residual sector reduces to:

$$
F_2(d)\neq 0
\quad\Longleftrightarrow\quad
(\mathrm{Re}\,F_2,\ \mathrm{Im}\,F_2)\neq (0,0),
$$

and it depends only on the **six** opposite-pair sums $S_k$.

### 6.2 Instantiation on exported coefficient-class language (R18 → P426)

`R18` already exports the six opposite-pair sums as explicit coefficient classes:

$$
\Sigma_{0,6},\ \Sigma_{1,7},\ \Sigma_{2,8},\ \Sigma_{3,9},\ \Sigma_{4,10},\ \Sigma_{5,11},
$$

where, for example,

$$
\Sigma_{0,6}
=
\left((3g4_{\psi0}v_{\psi0}^2+5g6_{\psi0}v_{\psi0}^4+2gY_0v_\phi^2+m2_{\psi0})-m_0^2\right)
+
\left((3g4_{\psi6}v_{\psi6}^2+5g6_{\psi6}v_{\psi6}^4+2gY_6v_\phi^2+m2_{\psi6})-m_0^2\right),
$$

with the certified floor $m_0^2 = 1.013551972358388$ (`R15`).

`P426` persists the exact reduced expressions for $\mathrm{Re}\,F_2$ and $\mathrm{Im}\,F_2$ in terms of these six
exported classes.

Concretely, identifying:

$$
S_0=\Sigma_{0,6},\ S_1=\Sigma_{1,7},\ S_2=\Sigma_{2,8},\ S_3=\Sigma_{3,9},\ S_4=\Sigma_{4,10},\ S_5=\Sigma_{5,11},
$$

the reduction becomes a direct 6-class linear combination:

$$
\mathrm{Re}\,F_2
=
\Sigma_{0,6}
\;+\;\frac{1}{2}\Sigma_{1,7}
\;-\;\frac{1}{2}\Sigma_{2,8}
\;-\;\Sigma_{3,9}
\;-\;\frac{1}{2}\Sigma_{4,10}
\;+\;\frac{1}{2}\Sigma_{5,11},
$$

$$
\mathrm{Im}\,F_2
=
\frac{\sqrt{3}}{2}\left(\Sigma_{1,7}+\Sigma_{2,8}-\Sigma_{4,10}-\Sigma_{5,11}\right).
$$

And the induced `pair1` diagonal `O(2)`-cut signature is:

$$
a_1-d_1=\frac{1}{6}\mathrm{Re}\,F_2,
\qquad
b_1=\frac{1}{12}\mathrm{Im}\,F_2.
$$

### 6.3 Two strict negative closures (what does *not* cut `O(2)` on `pair1`)

Two strict structural closures now eliminate common candidate stories for a diagonal/local `pair1` `O(2)` cut:

1. **FR parity sign mask cannot cut `O(2)` diagonally on `pair1`:** for the $n=12$ parity mask
   $b_k:=(-1)^k$, one has $F_2(b)=0$ (`N469`, audited in `P428`), hence no diagonal `pair1` `O(2)` cut can be supplied
   by a profile of the form $a+c\,(-1)^k$.
2. **+3 shift invariance kills the mode‑2 defect:** if a diagonal profile satisfies $d_{k+3}=d_k$ (indices mod $12$),
   then $F_2(d)=0$ (`N470`, audited in `P429`), hence it cannot cut `O(2)` on `pair1`.

## 7) The “most honest next move” after Release 6.4 (strict)

Release 6.4 turns a vague selector slogan into a single strict checkable question:

```text
Does the canonical FIN local diagonal residual sector have a nonzero mode-2 defect F2(d)?
```

On the current repo state this question is **not** strictly decidable from the exported canonical coefficient class
alone: `N472/P431` exhibit underdetermination (both outcomes remain compatible with current exports).
Additionally, `P432` confirms that no exported numeric instantiation of the opposite‑pair sums is currently present,
so even the reduced `P426` expression cannot yet be evaluated as a strict decision.

To proceed in strict core (without false pass), the repo must export at least one of:

1. **Strict-derived diagonal coefficient/value instantiation** sufficient to decide $F_2(d)\neq 0$, or
2. a strict theorem proving $F_2(d)=0$ (which would kill the diagonal-accelerator route on `pair1`), or
3. a different strict internal selector source (not a diagonal/local profile) that breaks translation invariance and
   cuts `O(2)` on `pair1` without introducing hidden slots.

Until such an ingredient is exported, no strict-core `pair1` `O(2)` cut, no strict-core theta export, and no
`QW-2191` discharge can be claimed.

## 8) Appendix: exported `pair1/pair2` mode-basis vectors on the 12-slot carrier (data)

For completeness, the repo exports explicit numeric vectors for the canonical 12-slot real Fourier mode scaffold
`(c1,s1,c2,s2)` as columns of the certified transport packet:

- `fundamental_action_reconstruction/generated/r11_symmetry_certified_declared_control_transport_packet_for_psi_block_route.json`

The values (in `psi0..psi11` order) are:

```json
{
  "c1": [0.408248290463863, 0.353553390593274, 0.204124145231932, 0.0, -0.204124145231931, -0.353553390593274, -0.408248290463863, -0.353553390593274, -0.204124145231932, 0.0, 0.204124145231932, 0.353553390593274],
  "s1": [0.0, 0.204124145231931, 0.353553390593274, 0.408248290463863, 0.353553390593274, 0.204124145231931, 0.0, -0.204124145231931, -0.353553390593274, -0.408248290463863, -0.353553390593274, -0.204124145231932],
  "c2": [0.408248290463863, 0.204124145231932, -0.204124145231931, -0.408248290463863, -0.204124145231932, 0.204124145231932, 0.408248290463863, 0.204124145231932, -0.204124145231931, -0.408248290463863, -0.204124145231931, 0.204124145231931],
  "s2": [0.0, 0.353553390593274, 0.353553390593274, 0.0, -0.353553390593274, -0.353553390593274, 0.0, 0.353553390593274, 0.353553390593274, 0.0, -0.353553390593274, -0.353553390593274]
}
```

---

## WERSJA POLSKA

## 0) Zmiana względem Release 6.3

Release 6.4 dodaje **nowe ścisłe domknięcia strukturalne** wokół przeszkody unikatowości `QW-2191`:

1. certyfikowany, translacyjnie niezmienniczy operator hosta jest **ściśle izotropowy** na
   `pair1 = span{c1,s1}` i nie może ciąć rodziny `O(2)` (`N465`, audit `P424`),
2. sektor diagonalny/lokalny tnie `O(2)` na `pair1` **wtedy i tylko wtedy**, gdy jego profil diagonalny ma niezerowy
   współczynnik Fouriera **trybu 2** `F2(d)` (`N466`, audit `P425`),
3. dla `n=12` defekt trybu 2 redukuje się do **6 sum par przeciwległych**
   `S_k := d_k + d_{k+6}` (dokładnie te 6 klas, które już eksportuje `R18`), więc ścisłe pytanie o cięcie `O(2)` na
   `pair1` redukuje się do jednego jawnego wyrażenia zespolonego (`N467`, artefakt `P426`).
4. jeśli `F2(d)≠0`, to oś/diagionalizacja na `pair1` jest kanoniczna i wyliczalna bez dodatkowych slotów:
   $\theta_*=\frac12\,\mathrm{atan2}(\mathrm{Im}\,F_2,\ \mathrm{Re}\,F_2)$ (`N468`, audit `P427`),
5. ścisła, FR‑pochodna maska parzystości $\mathbb{Z}_2$ $b_k:=(-1)^k$ ma **zerowy** defekt trybu 2:
   $F_2(b)=0$, więc sama nie może ciąć `O(2)` na `pair1` (`N469`, audit `P428`),
6. każdy profil diagonalny niezmienniczy względem przesunięcia o +3 ($d_{k+3}=d_k$ dla `n=12`) ma **zerowy** defekt
   trybu 2 $F_2(d)=0$, więc nie może ciąć `O(2)` na `pair1` (`N470`, audit `P429`).
7. dla `n=12`, jeżeli profil diagonalny jest niezmienniczy względem `Aut(Z_12)` (stały na 6 orbitach ilorazowych z
   `N455`), to defekt trybu 2 jest wymuszony jako **rzeczywisty** i redukuje się do jednego jawnego połączenia liniowego
   6 wartości orbit; więc pytanie o cięcie `O(2)` redukuje się do jednego rzeczywistego testu (`N471`, audit `P430`).
8. na obecnie wyeksportowanej kanonicznej klasie współczynników diagonalnych FIN (`R15`), `F2(d)` dla kanonicznego
   sektora `D_local_residual` jest **niedookreślone**: zarówno `F2(d)=0` jak i `F2(d)≠0` pozostają zgodne z eksportami
   (`N472`, audit `P431`).
9. repo nie eksportuje obecnie żadnej **instancji wartości liczbowych** 6 sum par przeciwległych (ani równoważnych
   wartości diagonalnych) wystarczającej do bezpośredniej ewaluacji `F2(d)` (`P432`).
10. defekt trybu 2 można odzyskać równoważnie z układu wpisów `R18` dla zadeklarowanego bloku `pair1`
   (`N473`, audit `P433`), a jawny harness ewaluacyjny jest przygotowany na moment, gdy pojawią się strict-derived
   wartości 6 sum par przeciwległych (`P434`).
11. przy stacjonarności stałej próżni (wyeksportowane kanoniczne EoM) i `vpsi_k ≠ 0`, kombinacja Yukawa znika z
    kanonicznego wpisu diagonalnego definiującego `d_k`, dając Yukawa‑wolną postać wyrażenia diagonalnego (`N474`,
    toy-audit `P435`).

To jest *redukcja*, nie rozładowanie: sprawia, że hasło “fizyczny akcelerator wyboru” staje się sprawdzalnym defektem,
ale nie dostarcza jeszcze strict-derived wartości współczynników diagonalnych. W szczególności `N472/P431` pokazują, że
kanoniczna klasa współczynników nie rozstrzyga `F2(d)` na obecnych eksportach, a `P432` potwierdza brak instancji
wartości liczbowych pozwalających rozstrzygnąć to przez bezpośrednią ewaluację.

## 1) Jednostronicowy stan ścisły (6.4)

### 1.1 Co jest teraz ostre strukturalnie

1. `QW-2191` pozostaje prawdziwe: sam kernel nie ustala kanonicznie wyboru bazy w zdegenerowanej podprzestrzeni 2D
   (ciągła rodzina `O(2)`).
2. Host tego nie rozwiązuje na `pair1`: certyfikowany operator hosta
   `A_host := K_total + m0^2 I` jest izotropowy na `pair1` (`N465`).
3. Każde ścisłe cięcie `O(2)` na `pair1` musi pochodzić z **nietranslacyjnie niezmienniczej** struktury, najnaturalniej
   z sektora diagonalnego/lokalnego z `F2(d) ≠ 0` (`N466`).
4. Dla kanonicznego sektora diagonalnego FIN na nośniku 12-slotowym, warunek `F2(d) ≠ 0` redukuje się do 6 klas z `R18`
   (`N467/P426`) – bez wprowadzania ukrytego “wyboru punktu” na pierścieniu.
5. Jeżeli takie cięcie diagonalne istnieje (`F2(d)≠0`), to oś/diagionalizacja na `pair1` jest kanoniczna i wylicza się
   wprost z `F2(d)` (`N468/P427`).
6. Dwa częste slogany o „łamaniu symetrii” są teraz zamknięte jako źródła diagonalnego cięcia `O(2)` na `pair1`:
   - maska parzystości FR $(-1)^k$ ma `F2=0` (`N469/P428`),
   - niezmienniczość względem przesunięcia o +3 ($d_{k+3}=d_k$) wymusza `F2=0` (`N470/P429`).
7. Jeżeli profil diagonalny jest dodatkowo niezmienniczy względem `Aut(Z_12)` (warunek ilorazowo-bezpieczny, bez ukrytych
   slotów), to `Im(F2)=0`, a kryterium diagonalnego cięcia `O(2)` na `pair1` redukuje się do jednego rzeczywistego testu
   na 6 wartościach orbit (`N471/P430`).
8. Na obecnie wyeksportowanej kanonicznej klasie współczynników diagonalnych FIN (`R15`), `F2(d)` dla kanonicznego
   sektora `D_local_residual` jest **niedookreślone** na poziomie eksportów, więc nie można jeszcze promować ścisłego
   diagonalnego świadka cięcia `O(2)` na `pair1` (`N472/P431`).

### 1.2 Co pozostaje brakujące (prawdziwy frontier)

Repo nadal **nie** eksportuje:

1. strict-derived instancji wartości/współczynników sektora diagonalnego wystarczającej do rozstrzygnięcia
   `F2(d) ≠ 0` (obecnie niedookreślone na eksportach: `N472/P431`),
2. ścisłego składnika theta/selektora kanonicznie tnącego `O(2)` (`T159`),
3. globalnego rozładowania `QW-2191` ani domknięcia ToE.

Release 6.4 jest więc krokiem redukcyjnym: zastępuje retorykę selektora jednym ścisłym defektem `F2(d)`.

## 2) Ścisły kernel

Aktualny ścisły kernel roboczy:

$$
K_{\mathrm{sg}}(d)
=
\frac{\cos(\omega d+\phi)}{1+\beta d^{\eta}},
\qquad
(\omega,\phi,\beta,\eta)=(0.18575,\ 0.16250,\ 1.0,\ 1.8).
$$

## 3) Nowy wynik ścisły A: izotropia hosta na `pair1` (N465 + P424)

W bazie `(c1,s1)` audit `P424` daje:

$$
\left[A_{\mathrm{host}}\right]_{\{c_1,s_1\}}
=
\begin{pmatrix}
1.9197380969174072 & 5.144911710811878\times 10^{-17}\\
5.144911710811878\times 10^{-17} & 1.9197380969174067
\end{pmatrix},
$$

co jest zgodne z izotropią w tolerancji $10^{-12}$.

## 4) Nowy wynik ścisły B: kryterium cięcia `O(2)` przez sektor diagonalny (N466)

Dla $D=\mathrm{diag}(d_0,\ldots,d_{n-1})$:

$$
F_2(d) := \sum_{i=0}^{n-1} d_i\, e^{i\frac{4\pi i}{n}},
$$

a sektor diagonalny tnie `O(2)` na `pair1` wtedy i tylko wtedy, gdy $F_2(d)\neq 0$.

## 5) Nowy wynik ścisły C: redukcja defektu trybu 2 do 6 klas (N467 + P426)

Dla `n=12`:

$$
F_2(d)=\sum_{k=0}^{5} S_k e^{i\frac{\pi k}{3}},
\qquad
S_k := d_k + d_{k+6}.
$$

W szczególności:

$$
\mathrm{Re}\,F_2
=
(S_0 - S_3) + \frac{1}{2}(S_1 - S_2 - S_4 + S_5),
$$

$$
\mathrm{Im}\,F_2
=
\frac{\sqrt{3}}{2}(S_1 + S_2 - S_4 - S_5).
$$

`R18` eksportuje dokładnie te 6 sum jako klasy `\Sigma_{k,k+6}`, a `P426` zapisuje jawne wyrażenia dla
`Re(F2)` i `Im(F2)` wprost w tym języku.

Konkretnie (dla identyfikacji `S_0=\Sigma_{0,6}`, ..., `S_5=\Sigma_{5,11}`):

$$
\mathrm{Re}\,F_2
=
\Sigma_{0,6}
\;+\;\frac{1}{2}\Sigma_{1,7}
\;-\;\frac{1}{2}\Sigma_{2,8}
\;-\;\Sigma_{3,9}
\;-\;\frac{1}{2}\Sigma_{4,10}
\;+\;\frac{1}{2}\Sigma_{5,11},
$$

$$
\mathrm{Im}\,F_2
=
\frac{\sqrt{3}}{2}\left(\Sigma_{1,7}+\Sigma_{2,8}-\Sigma_{4,10}-\Sigma_{5,11}\right).
$$

Oraz:

$$
a_1-d_1=\frac{1}{6}\mathrm{Re}\,F_2,
\qquad
b_1=\frac{1}{12}\mathrm{Im}\,F_2.
$$

### 5.1 Dwa ścisłe domknięcia negatywne (co *nie* tnie `O(2)` na `pair1`)

Dwa ścisłe fakty strukturalne eliminują częste kandydackie historie o diagonalnym/lokalnym cięciu `O(2)` na `pair1`:

1. **Maska parzystości FR nie tnie `O(2)` diagonalnie na `pair1`:** dla maski $n=12$
   $b_k:=(-1)^k$ zachodzi $F_2(b)=0$ (`N469`, audit `P428`), więc profil postaci $a+c\,(-1)^k$ nie może sam dostarczyć
   diagonalnego cięcia `O(2)` na `pair1`.
2. **Niezmienniczość względem przesunięcia o +3 zabija defekt trybu 2:** jeżeli profil diagonalny spełnia
   $d_{k+3}=d_k$ (indeksy mod $12$), to $F_2(d)=0$ (`N470`, audit `P429`), więc nie może ciąć `O(2)` na `pair1`.

## 6) Najuczciwszy następny ruch (strict)

Release 6.4 sprowadza pytanie o “akcelerator wyboru” do jednego ścisłego testu:

```text
czy kanoniczny sektor diagonalny FIN ma niezerowy defekt trybu 2: F2(d) ≠ 0 ?
```

Na obecnym stanie repo test ten nie jest rozstrzygalny w ścisłym rdzeniu z samej wyeksportowanej kanonicznej klasy
współczynników: `N472/P431` pokazują niedookreślenie (oba wyniki są kompatybilne z eksportami).
Dodatkowo `P432` potwierdza brak jakiejkolwiek wyeksportowanej instancji wartości liczbowych 6 sum par przeciwległych,
więc nawet zredukowane wyrażenie `P426` nie może zostać dziś policzone jako ścisła decyzja.

Żeby iść dalej w ścisłym rdzeniu (bez fałszywego PASS), repo musi więc wyeksportować strict-derived wartości/relacje
wystarczające do rozstrzygnięcia tego testu (albo dowód `F2(d)=0`).

## 7) Aneks: wyeksportowane wektory bazy trybów na nośniku 12-slotowym (dane)

Repo eksportuje jawne wartości wektorów `(c1,s1,c2,s2)` (kolejność `psi0..psi11`) w:

- `fundamental_action_reconstruction/generated/r11_symmetry_certified_declared_control_transport_packet_for_psi_block_route.json`

```json
{
  "c1": [0.408248290463863, 0.353553390593274, 0.204124145231932, 0.0, -0.204124145231931, -0.353553390593274, -0.408248290463863, -0.353553390593274, -0.204124145231932, 0.0, 0.204124145231932, 0.353553390593274],
  "s1": [0.0, 0.204124145231931, 0.353553390593274, 0.408248290463863, 0.353553390593274, 0.204124145231931, 0.0, -0.204124145231931, -0.353553390593274, -0.408248290463863, -0.353553390593274, -0.204124145231932],
  "c2": [0.408248290463863, 0.204124145231932, -0.204124145231931, -0.408248290463863, -0.204124145231932, 0.204124145231932, 0.408248290463863, 0.204124145231932, -0.204124145231931, -0.408248290463863, -0.204124145231931, 0.204124145231931],
  "s2": [0.0, 0.353553390593274, 0.353553390593274, 0.0, -0.353553390593274, -0.353553390593274, 0.0, 0.353553390593274, 0.353553390593274, 0.0, -0.353553390593274, -0.353553390593274]
}
```
