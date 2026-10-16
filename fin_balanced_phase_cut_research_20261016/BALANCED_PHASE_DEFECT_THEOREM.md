# New paid cut class: bounded-depth defects of a balanced phase

Date: 2026-10-16. Exact finite analytical estimate under actual local native
ellipticity and phase symmetry. Not an all-cut bound or constructed kinetics.

## 1. Phase-centering removes spurious bulk terms

Keep s=A_U u, |s_x|<=N W_x, exact s(A)=s(B)=0, and let

    S=(A\T_minus) union T_plus,
    T_minus subset A, T_plus subset B.

Then **s(S)=s(T_plus)-s(T_minus)** exactly. There is no unsigned whole-phase
interface term in this identity. Every internal native edge has

    Rw^-1 <= W_y/W_x <= Rw,
    Rw=p_max/a0>=1.

This follows by applying the lower/upper native rate bounds to both
orientations and detailed balance. These are reference-weight inequalities,
not assumptions on a diagnostic generator.

## 2. Oriented defect-depth hypothesis and conclusion

Suppose a fixed integer L has the following paid properties:

- Every x in T_minus has an upward H path v->...->x of length <=L, whose
  first vertex v is in S and all remaining vertices are in T_minus.
- Every x in T_plus has an upward H path x->...->v->y of length <=L, whose
  vertices except y lie in T_plus and whose last vertex y lies in S^c.

Paths may be chosen deterministically for the proof. Every actual upward
count vertex has at most 36 successors and 36 predecessors. Charge a missing
vertex to its path's entrance boundary edge, and an added vertex to its
path's exit boundary edge. For each edge, the number of distinct charged
vertices is at most `1+36+...+36^(L-1)<36^L`, separately for each kind.
Neighbor-weight comparison gives W_x<=Rw^L W_tail(boundary edge). Hence

    W(T_minus)+W(T_plus)<=2(36Rw)^L B(S),
    [s(S)]_+/(N B(S))<=2(36Rw)^L.                    (1)

Zero denominators are included: the path assumptions force empty defects
or zero total charged weight when B=0. The constant is very conservative
but independent of N. It is not an assertion that every cut has bounded
depth; that assertion would be a new, generally unsupported hypothesis.

## 3. An actual automatic subcase

Suppose the defects comprise at most J original count states in total, all
outside the true departure E_N. The earlier phase geometry supplies a
within-phase H predecessor and successor at every such state, for sufficiently
large N. Repeatedly follow predecessors within T_minus until leaving it;
strict decrease of m prevents repetition. At most J steps reach A\T_minus.
Reverse that path. Similarly, repeated successors from T_plus leave it
within J steps, in B\T_plus. Thus the depth hypothesis holds with L<=J.

If at most J **G_+ orbits** are changed, their total original size is at
most 12J, giving L<=12J. No assumption of free action or equal orbit sizes
is needed. This proves a uniform score bound for every bounded-orbit
interior perturbation of a balanced phase. The actual cuts used to disprove
the stronger shielding lemma are among these benign cuts.

## 4. Remaining dependency, not a new axiom

A diverging original cut score cannot be established merely by a bounded
interior orbit perturbation of a whole phase. It must escape this depth
class, involve growing phase-defect depth, or require genuine departure-
layer treatment. This is a necessary restriction on proposed adverse cuts,
not a classification of every possible one.

For general cuts the algebraic phase-centering identity still holds. What
is unpaid is a **source-sensitive**, signed estimate for long-depth defects
and far soft tails, with actual boundary permissions. The rough count-depth
constant in (1) grows exponentially with L and is not a uniform all-N bound
when L grows. Do not insert bounded depth, expansion or a lower diagnostic
activity as a microscopic law to conceal this remaining estimate.
