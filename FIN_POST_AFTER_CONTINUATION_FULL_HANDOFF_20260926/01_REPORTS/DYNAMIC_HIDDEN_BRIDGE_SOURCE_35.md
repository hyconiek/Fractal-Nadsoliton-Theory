# DYNAMIC-HIDDEN-BRIDGE-SOURCE-35
## Exact uniqueness in the minimal A–D_u grammar, and exact nonuniqueness beyond it

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`,
visible HEAD `fe14a6f4e436815635df54429102f23c22862296`.

Status:
- **EXACT algebraic theorem** inside a declared minimal operator grammar;
- **EXACT no-go** for uniqueness under the broader endpoint-local D12 class;
- physical provenance remains open.

This continues `RELATIONAL-CLOCK-HIDDEN-BRIDGE-34`.

---

## 1. Existing typed ingredient

The repository already contains the conditional state-to-multiplication object

\[
u\mapsto D_u=\operatorname{diag}(u),
\]

covariant under D12:

\[
D_{Pu}=P D_u P^T.
\]

The strict operator \(A\) is real symmetric, D12-invariant, and

\[
A\mathbf1=0.
\]

The question is whether these already admitted ingredients force a
hidden-to-generator bridge.

---

## 2. Minimal operator grammar

Consider the smallest first-order matrix grammar that is

- linear in \(u\),
- linear in one occurrence of \(A\),
- built from ordinary matrix multiplication,
- allowed one diagonal correction formed from \(Au\).

The general element is

\[
J(u)
=
a\,D_uA
+b\,AD_u
+c\,\operatorname{diag}(Au).
\]

No shell projectors, entrywise functions of \(A\), extra kernels or fitted
coefficients are admitted.

---

## 3. Symmetry forces \(a=b\)

Because \(A=A^T\) and \(D_u=D_u^T\),

\[
(D_uA)^T=AD_u.
\]

Demanding

\[
J(u)^T=J(u)
\]

for every \(u\) therefore gives

\[
a=b.
\]

Write this common value as \(\alpha\).

---

## 4. Conservation forces the diagonal correction

Using

\[
A\mathbf1=0,\qquad D_u\mathbf1=u,
\]

we have

\[
J(u)\mathbf1
=
\alpha Au+cAu
=
(\alpha+c)Au.
\]

For nonconstant hidden modes \(u\in H_4\), \(Au\neq0\).

Thus preservation of the constant/zero mode,

\[
J(u)\mathbf1=0,
\]

forces

\[
c=-\alpha.
\]

Therefore:

\[
\boxed{
J(u)
=
\alpha\bigl(
D_uA+AD_u-\operatorname{diag}(Au)
\bigr)
}
\]

is the **unique bridge up to one scalar** inside this minimal grammar.

The normalization used in report 34 corresponds to \(\alpha=1/2\).

---

## 5. D12 covariance is automatic

Let \(P\) be any D12 permutation.

Because

\[
PA=AP,\qquad
D_{Pu}=P D_u P^T,
\]

the unique minimal-grammar bridge satisfies

\[
J(Pu)=P J(u)P^T.
\]

Thus no additional covariance tuning is required.

---

## 6. Dirichlet interpretation

For \(\alpha=1/2\), the off-diagonal entries are

\[
J_{ij}
=
\frac12 A_{ij}(u_i+u_j),\qquad i\ne j.
\]

Since \(A_{ij}=-w_{ij}\), this corresponds exactly to

\[
\delta w_{ij}
=
\frac12 w_{ij}(u_i+u_j).
\]

So the minimal algebraic bridge equals the first variation of a Dirichlet
network whose edge conductance is modulated by the endpoint-average field.

This gives the bridge a clean variational interpretation.

However, the statement that the physical hidden field must modulate the
Dirichlet density through its endpoint average is still an additional
constitutive premise.

---

## 7. Why this is not yet a FIN source theorem

Drop only one restriction: allow the bridge to resolve the six unoriented
cyclic edge shells.

For shell \(d=1,\dots,6\), let \(w^{(d)}_{ij}\) be the base weight restricted to
that shell.

Then every rule

\[
\delta w_{ij}
=
c_d\,(u_i+u_j)
\]

is

- linear in \(u\),
- endpoint-local,
- symmetric,
- D12-equivariant,
- Laplacian/conservation preserving.

Treating these as linear maps from the complete hidden space \(H_4\) to
symmetric zero-row-sum operators, the six shell maps are linearly independent.

Therefore the broader admissible space has dimension exactly six:

\[
(c_1,\dots,c_6)\in\mathbb R^6.
\]

Examples such as "modulate shell 1 only" and "modulate shell 2 only" satisfy
the same broad symmetry/conservation/locality requirements and are inequivalent.

Hence:

\[
\boxed{
\text{D12 + symmetry + conservation + endpoint locality do not select }J.
}
\]

The one-parameter bridge becomes unique only after an extra restriction is
declared.

---

## 8. Exact location of the missing law

At least one of the following must be sourced:

### Route A — minimal algebra grammar
Only \(A\), \(D_u\), ordinary composition and the required diagonal
conservation correction are admissible at first order.

### Route B — shell-blind fractional response
All shells obey the same fractional law

\[
\delta w_{ij}/w_{ij}
\propto u_i+u_j.
\]

### Route C — a variational continuum/refinement principle
A refinement or action law derives endpoint-average modulation as the unique
discretization compatible with the declared limit.

Without one of these, the bridge remains one member of a six-dimensional
constitutive family.

---

## 9. Relation to existing ST233 state-to-multiplication result

ST233 establishes that once a state/field is supplied,

\[
u\mapsto D_u
\]

is a genuine typed second operator and need not commute with \(A\).

The present result adds:

- if one insists on the minimal one-\(A\), one-\(D_u\) symmetric conservative
  grammar, the induced Laplacian variation is unique up to scale;
- ST233 itself does not justify that grammar;
- therefore localization/state information is available, but its dynamical
  coupling remains a separate constitutive choice.

This is exactly consistent with the earlier source-boundary discipline.

---

## 10. Consequence for emergent physics

The research chain now distinguishes:

1. **hidden state exists mathematically** — four discarded modes;
2. **hidden state is observable in coarse dynamics** — quadratic detector theorem;
3. **a natural hidden-to-generator bridge exists** — report 34;
4. **that bridge is unique in a minimal algebraic grammar** — this report;
5. **the grammar itself is not yet sourced** — six-shell no-go.

So the open problem is narrower than "invent dynamics".

It is now:

> Why should nature/FIN use the total operator \(A\) as an indivisible coupling
> object, rather than allowing shell-resolved constitutive responses?

That is a concrete source question.

---

## 11. Next research atom

### REFINEMENT-FORCES-BRIDGE-36

Test whether the already admitted subdivision/refinement laws can eliminate the
six-shell freedom.

For each shell-resolved coefficient vector \(c=(c_1,\ldots,c_6)\):

1. refine an edge/network according to the accepted static/dynamic refinement
   rules;
2. require coarse-graining to commute with first-order hidden modulation;
3. solve the resulting linear constraints on \(c_d\).

Acceptance:

- the constraint space collapses to
  \[
  c_d\propto w_d,
  \]
  deriving the fractional bridge up to scale; or
- dimension \(>1\) survives, proving refinement alone cannot source the bridge.

This is the next high-value physical-source test.
