# FIN ST3822–ST3896 — Five-Round Hausdorff-Moment Obstruction

## Question

Can the six positive strict conductances be the beginning of a completely
monotone relaxation law

\[
K(d)=\int_0^1x^d\,d\mu(x),
\qquad \mu\ge0?

\]

Such a representation would provide a positive mixture of decaying
exponentials and a distinguished smooth local continuation.

## Round 1: deceptive finite-difference pass

All available alternating finite differences of
\(K(1),\ldots,K(6)\) are positive. The elementary complete-monotonicity screen
therefore passes. This condition is necessary but insufficient for the
truncated moment problem.

## Round 2: Hankel obstruction

A positive representing measure would make the shifted moment matrices

\[
H_1=(K(i+j+1))_{i,j=0}^2,
\qquad
H_2=(K(i+j+2))_{i,j=0}^2

\]

positive semidefinite, because they are Gram matrices of powers of \(x\).

The strict matrices instead have

\[
\lambda_{\min}(H_1)=-0.000565002186446,

\]

\[
\lambda_{\min}(H_2)=-0.001100659961814.

\]

Their determinants are negative. Therefore:

\[
\boxed{
\text{no positive measure on }[0,1]
\text{ reproduces all six strict weights}.}

\]

Equivalently, there is no exact positive mixture of real decaying exponentials
with those moments.

## Round 3: robustness

If every moment is perturbed by at most \(\varepsilon\), the corresponding
three-by-three matrix perturbation has norm at most \(3\varepsilon\). Repairing
both negative eigenvalues therefore requires at least

\[
\boxed{\varepsilon\ge0.000366886653938.}

\]

This is approximately 3.31 percent of \(K(6)\). Ordinary floating-point or
decimal rounding is far too small to remove the obstruction. An approximate
positive mixture requires a material deformation and a chosen loss function.

## Round 4: legacy comparison

The canonical legacy sequence contains negative values and immediately fails
positive moment representation. Taking absolute values does not repair it:
both tested Hankel matrices remain strongly indefinite. Signed or complex
exponential representations remain possible but lose the positive-relaxation
interpretation and are nonunique.

## Round 5: verdict

\[
\boxed{3/10\ \text{mathematical rows},\qquad0/10\ \text{physical rows}.}
\]

The strict finite kernel is positive as a graph-conductance sequence on its
declared six distances, but it is not a truncated Hausdorff moment sequence.
Consequently a positive Stieltjes/Bernstein or exponential-mixture tail cannot
be the missing canonical continuum without changing the frozen weights.

The next viable route must classify positive local continuations under exact
finite matching by a principle other than complete monotonicity, or explicitly
declare an approximate deformation and its loss.

No continuum selection, physical units, apparatus evidence, Standard Model,
gravity, \(L_{\rm total}\), or Theory of Everything closure follows.
