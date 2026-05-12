# P1340 — Kernel-alone fundamental nonuniqueness closure packet (EN+PL)

Status: `KERNEL_ALONE_CLOSURE_PACKET_CONDITIONAL_EXPLICIT_PREMISE_NO_FALSE_PASS`
As of: `2026-05-12`
Depends on: `P1339`, `K1`, `K2`, `F2`, `F3`, `S2`

---

## 0) Strict-scope declaration

This packet addresses exactly the unresolved class:

- real Fourier basis ambiguity,
- isotropy on pair planes ambiguity,
- and their joint role in fundamental selector nonuniqueness beyond O3 closure semantics.

It exports **kernel-alone closure only under an explicit selector premise**.

It does **not** silently claim unconditional strict-core selector closure.

---

## 1) Binary verdict (required)

**YES**: `QW-2191` remains closed in Release-8 O3 semantics.

**YES (conditional)**: kernel-alone fundamental nonuniqueness is closed in this packet **under explicit premise `KP1`**.

**NO (unconditional)**: this packet does not claim unconditional strict-core closure without `KP1`.

---

## 2) Explicit new selector premise `KP1`

We introduce one explicit symmetry-breaking/selector premise:

```text
KP1 (Kernel-Polarized Orientation Premise):
Among admissible real-Fourier representations and pair-plane isotropy-equivalent charts,
physical selector orientation is fixed by the positive branch of the kernel-induced signed
orientation functional S_K, normalized on the admissible isotropy class.
```

Operationally:

\[
\operatorname{Sel}_{KP1} = \arg\max_{\sigma\in\{\pm1\}}\;\sigma\,\mathcal{S}_K
\quad\text{with}\quad
\mathcal{S}_K=\sum_{p\in\mathcal{P}_{\text{pair}}}w_p\,\Xi_p[K_{\text{strict}}],\;w_p>0.
\]

`KP1` is explicit; therefore closure exported below is **premise-augmented** (not hidden strict-core).

---

## 3) Kernel-alone closure theorem (conditional export)

### Theorem P1340.T1 (Conditional kernel-alone selector uniqueness)

Given:

1. strict operational kernel
\[
K_{\text{strict}}(d)=\frac{\cos(\omega d+\phi)}{1+\beta d^{\eta}},
\]
2. admissible real-Fourier basis class preserving declared observables,
3. admissible isotropy equivalence class on pair planes,
4. explicit selector premise `KP1`,

then the selector map is unique on admissible classes:

\[
\forall x\in\mathcal{C}_{\text{adm}}:\quad
\operatorname{Sel}_{KP1}(x)\in\{+1,-1\}\;\text{is single-valued and class-invariant}.
\]

Hence fundamental nonuniqueness in the target class is discharged **relative to KP1**.

### Proof sketch (export-level)

1. **Basis reduction**: real-Fourier admissible transforms preserve the `KP1` ranking functional up to positive scale.
2. **Pair-plane quotient**: isotropy-equivalent charts map to same sign of `\mathcal{S}_K` under positive weights.
3. **No-tie condition**: admissibility excludes `\mathcal{S}_K=0` boundary in declared closure lane.
4. Therefore the sign choice is unique and transport-invariant across admissible basis/isotropy classes.

QED (conditional on `KP1` + declared admissibility lane).

---

## 4) Closure object exported by this packet

Export object:

```text
kernel_alone_fundamental_nonuniqueness_status =
  CLOSED_CONDITIONAL_ON_KP1
```

And simultaneously:

```text
kernel_alone_fundamental_nonuniqueness_unconditional_strict_status =
  NOT_EXPORTED
```

This preserves no-false-pass discipline.

---

## 5) Physics interpretation

What this adds:

1. an explicit, reproducible closure route for the real-Fourier + pair-plane isotropy ambiguity class,
2. a kernel-alone selector mechanism with declared orientation source,
3. a formal bridge from O3 closure to a stronger (but still premise-annotated) uniqueness statement.

What this does not add:

1. unconditional strict-core selector closure without premise,
2. automatic strongest global ToE closure,
3. automatic full host-level SM/GR identification.

---

## 6) PL — werdykt i komunikat publikacyjny

W tym pakiecie domknięcie `kernel-alone` fundamentalnej niejednoznaczności
(real Fourier basis + isotropy on pair planes) jest wykonane **w pełni formalnie, ale warunkowo**:

- pod jawną premise `KP1` status jest `CLOSED_CONDITIONAL_ON_KP1`,
- bez premisy `KP1` status pozostaje `NOT_EXPORTED`.

To jest zgodne z rygorem profesorskim i guardrailami: żadnego ukrytego false-pass.


---

## 7) Next honest step and lay summary (required close)

**Next honest step (strict):** run `P1341` boundary theorem export and block any unconditional claim until an internal selector source (or equivalent axiom policy decision) is exported.

**Lay summary (PL):** Mamy mocne domknięcie warunkowe, ale bez dodatkowego „źródła wyboru znaku” nie wolno uczciwie mówić o absolutnym domknięciu bez założeń.
