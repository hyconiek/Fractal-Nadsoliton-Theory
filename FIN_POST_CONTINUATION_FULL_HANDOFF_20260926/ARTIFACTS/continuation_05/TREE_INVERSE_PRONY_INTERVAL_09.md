# TREE-INVERSE-CONDITIONING-09 — explicit two-pole interval certificate

Status: **PROOF-GRADE SUFFICIENT INTERVAL CERTIFICATE**.

For a scalar projection of a passive response with two visible poles,
  m_n=a1 lambda1^n+a2 lambda2^n,
let
  D=m0 m2-m1^2=a1 a2 (lambda1-lambda2)^2.
The Prony symmetric coefficients are
  s1=(m0 m3-m1 m2)/D,
  s2=(m1 m3-m2^2)/D.

Given boxes |mhat_j-m_j|<=eps_j, define
  deltaD <= |mhat2|eps0+|mhat0|eps2+2|mhat1|eps1+eps0 eps2+eps1^2,
and analogous product-rule boxes for the two numerators.  If the lower bounds
of D, N1 and N2 are positive, interval division gives certified intervals for
s1,s2.  A sufficient separated-pole certificate is
  s1_lower^2 - 4 s2_upper > 0.

`prony_interval_certificate.py` evaluates this without fitting.  In the fixture
lambda=1+-delta, a1=a2=1, with independent absolute moment boxes equal to
1e-6 max(1,|m_n|), the conservative certificate requires approximately
  delta > 0.0669391,
i.e. pole separation >0.1338782.  The true inverse still exists below this;
the statement is only that this error box no longer certifies two separated
poles.

This quantifies the earlier conditioning identity and remains separate from the
static short-edge topology margin and the matrix-residue visibility margin.
