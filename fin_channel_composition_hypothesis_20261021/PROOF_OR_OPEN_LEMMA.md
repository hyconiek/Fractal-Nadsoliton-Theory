# Classification, adverse controls and actual finite-parent consequence

Date: 2026-10-21. Conditional analytical proofs, not formalization or
independent review. All premises are frozen in the authorization/contract.

## CC-001. Complete response classification

Put ell_ab(h)=log(q R_ab(h)). Positivity and channel composition make ell
additive; continuity at zero makes it a real linear functional. Write
ell_01(h)=sum_j c_j h_j. Permutations fixing 0 and 1 force all c_j, j>=2,
to be equal. Covariance under exchanging 0 and 1 and detailed balance give

    c_1-c_0=1.

Shift invariance gives sum c_j=0. For q>=3 this leaves exactly one real
parameter s, with

    R_ab^s(h)=(1/q)exp[(h_b-h_a)/2+s(h_a+h_b-2 mean(h))].       (1)

Conversely every fixed real s satisfies every premise. The coefficient
representation is unique for q>=3. For q=2 the symmetric expression
h_a+h_b-2 mean(h) is identically zero, and the unique rule is the square-root
rule. These are continuous off-diagonal rates, not transition probabilities.

The symmetric edge activity in (1) is

    a_ab^s(h)=sqrt(R_ab R_ba)=(1/q)exp[s(h_a+h_b-2 mean(h))].  (2)

Thus the new principle genuinely restricts the formerly arbitrary activity
functions, but it DOES NOT determine a unique activity for q12. Choosing
s=0, s=1/2 or any other value would be another choice, not a conclusion.
In particular no claim of a naturally selected HB or SQRT rule is paid.

An endpoint-only axiom would force s=0 for q>=3: varying any other h_j
changes the exponent in (1) by -2s h_j/q. That axiom is NOT part of this
model and is not adopted as a repair. Continuity, full rule covariance and
the neutral reference also cannot be silently deleted from the theorem.

## CC-002. Native HB and the mandatory thinning control

Native HB has R_ab(h)=exp(h_b)/Z(h), Z=sum exp(h_j). Its composition defect
is independent of a,b:

    log[R_ab(u+v)R_ab(0)/(R_ab(u)R_ab(v))]
      =log Z(u)+log Z(v)-log q-log Z(u+v).                    (3)

For u=v nonconstant this is strictly negative by strict Cauchy–Schwarz:
(sum exp(u))^2<q sum exp(2u). Equivalently, along h=t v,

    d^2/dt^2 log R_ab(t v)=-Var_softmax(tv)(v)<0.

In (1) that derivative is zero. Hence HB fails the response-rule axiom in
the frozen neutral clock. This failure is exact for any nonconstant actual
A7 conditional field; no rounded spectral or critical-gain value is needed.

For the prior distance-three-thinned HB family, the extra log factor is
-log(1+g^2) on those edges. Along h=g v its second derivative at zero adds
-2 to the HB curvature. It cannot restore composition. On other edges the
original strict HB defect remains. Both native and thinned HB remain valid
old same-Gibbs models; rejection here concerns membership in the NEW
postulated class, not their mathematical validity or nature.

## CC-003. The clock trap at N2 and its exact N3 resolution

For fixed h, destination odds under (1) are

    R_ab/R_ac=exp[(1/2+s)(h_b-h_c)].                         (4)

To match HB odds at two distinct allowed destination logits requires s=1/2.
At this value

    R_ab^(1/2)(h)=exp(h_b-mean(h))/q
                =[Z(h)exp(-mean(h))/q] R_ab^HB(h).          (5)

For a single supplied categorical field, this is only a common clock change.
In actual q12/N2 with a circulant row-zero A7, every one-background field
is a rotation of the same vector. Its Z and mean are independent of the
background label. Therefore s=1/2 DOES match the entire N2 HB generator
modulo one gain-dependent global factor. Do not use N2 alone to claim a
clock-invariant discrimination or to choose s=1/2 as a sourced value.

Now consider the actual labelled N3 state x=(0,0,1). A slot carrying 1
has field H=(2g/3)A7 e_0. A slot carrying 0 has field

    M=(g/3)A7(e_0+e_1)=(H+shift(H))/2.

Both means are zero. Both rotations of H have the same Z. The columns
A7 e_0 and A7 e_1 differ: otherwise their nonzero parity eigencomponent
would be equal, contrary to Lambda6>0 and opposite parity of 0 and 1.
Strict convexity of exp now gives

    Z(M)<Z(H), for every g>0.                               (6)

There is also a direct square identity, not just a convexity citation:

    Z(H)-Z(M)=(1/2)sum_j [exp(H_j/2)-exp(shift(H)_j/2)]^2>0.

To force s=1/2, the slot carrying 1 has two different allowed destination
logits. Indeed, if the reflection-symmetric row of A7 were constant off
label 1, then equality A_01=A_0,11 would force it to be constant everywhere,
contradicting Lambda6>0. Equation (4) therefore applies at this very state.

For s=1/2 the clock multipliers for jumps of the 1 slot and of a 0 slot
are respectively Z(H)/q and Z(M)/q. They are different by (6). Therefore
NO s in the classified family can match actual N3 HB even by one clock
factor depending on the FULL current state. A fortiori no single global
clock suffices. This is a finite labelled-law result; no MT premise or
large-N approximation enters it.

There is a bounded count-observable consequence. At n=(2,1,0,...), condition
on the next count change and ask whether n_1 decreases. Let a>0 be the HB
exit rate of the 1 slot and b>0 the combined HB exit rate of the two 0
slots. HB gives a/(a+b). The s=1/2 model gives

    Z(H)a/[Z(H)a+Z(M)b] > a/(a+b).                          (7)

The probability gap is exactly
`ab[Z(H)-Z(M)]/[(a+b)(Z(H)a+Z(M)b)]`, so its strict sign follows from
the displayed square and the positive legal exit rates.

This is an ideal conditional model comparison. It is not a measured
probability, a late-iid critical theorem or a claim about unspecified partial
readouts. Null refresh attempts are not observed count changes and are not
included as extra events. Other values of s already fail the destination-odds test.

For a separate exact arithmetic control, take A=q I-11^T and
g=3 log(2)/q at N3. Gibbs weights are proportional to 2^(equal-pair count).
The two conditional unnormalized rows have sums q+3 and q+2, with the same
product. At q12 equation (7) reads

    HB: 49/139,   s=1/2: 7/19,   difference: 42/2641>0.

This exact number belongs to the stated control, NOT to versioned A7.
Its harmless global algebraic clock cancels from the conditional event.

## CC-004. Actual Gibbs/count compatibility and clocks

Keep the actual FIN LOO field h=(g/N)A7(n-e_a). The background is identical
on the reverse edge; constant diagonal makes the labelled Gibbs potential
difference h_b-h_a. Equation (1) therefore has exact labelled detailed
balance. Slot covariance gives the quotient n_a R_ab(h); its reverse has
the count multiplicity n_b+1. Consequently it has the original count Gibbs
law with multinomial multiplicity, not a uniform count law.

At bounded |g|<=G and fixed s, |h_j|<=G M_A and mean(h)=0. Hence

    total exit <= N (q-1)/q exp[(1+2|s|)G M_A].              (8)

This supplies a finite N-independent cap for bounded gain. It does NOT
give a per-slot unit attempt clock. In fact no member can have a uniform
finite exit cap over all h: each ell_ab is linear and detailed balance
forces at least one to be nonzero; along one of the two rays h=+t v or
-t v that rate diverges. A unit-cap implementation valid for arbitrary
fields is incompatible with the frozen channel rule. Normalizing the
rates afterward changes the tested hypothesis.

## CC-005. What leading collective data cannot select

This is an exact mean-field coefficient calculation, not a new path theorem.
Use p_m=(1+m sigma)/q and h=x sigma, x=g Lambda6 m/q. For a parity-changing
transition h_a+h_b=0, so its rate is e^(-sigma_a x)/q for EVERY s.
The scalar mean-field parity drift is consequently

    sinh(x)-m cosh(x),                                     (9)

independent of s. At uniformity its per-copy parity jump second moment is
2. At g_c, equation (9) has cubic term -m^3/3. Same-parity exchanges do
depend on s; they cannot be inferred from this scalar drift.

More generally, for arbitrary zero-sum first-order perturbations w and h
at uniformity, differentiating inflow minus outflow in (1) gives

    B_linear=-w+h/q,                                      (10)

again independent of s. Thus a leading linear-response or parity-line test
does not fix the remaining activity. This conclusion must not be inflated
to equality of finite-N histories or an unproved critical process theorem
for every s. The old SQRT theorem is reusable only at s=0, with its explicit
factor-q clock conversion and its actual initialization assumptions.

## Terminal scope

The authorized second principle gives a real conditional reduction of the
activity freedom, but leaves a free parameter and fails as a provider of
the original N3 HB law even modulo state clocks. No s is selected or fitted.
Its independent physical justification, and the supplied A7/Gibbs/carrier/
clock premises, remain unpaid. Natural F01 remains OPEN_SOURCE. Return to
architecture rather than append an endpoint, cap or clock axiom as repair.
