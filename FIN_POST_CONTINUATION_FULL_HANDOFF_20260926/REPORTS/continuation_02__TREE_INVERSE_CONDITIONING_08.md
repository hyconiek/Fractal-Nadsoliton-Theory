# TREE-INVERSE-CONDITIONING-08

Status: PROOF-GRADE CONDITIONAL MATHEMATICS for the declared passive pole/residue response.

Let

Lambda(z) = L_BB - sum_r R_r/(z+lambda_r),   R_r >= 0, lambda_r>0.

For any boundary projection u define a_r=u^T R_r u >=0 and moments
m_n=sum_r a_r lambda_r^n.  These are the scalar high-frequency moments of the
response.

## Exact Hankel/Vandermonde conditioning identity

For the 2x2 Hankel minor,

D_2 = m_0 m_2 - m_1^2
    = sum_{r<s} a_r a_s (lambda_r-lambda_s)^2.

More generally, for H_q=[m_{i+j}]_(i,j=0..q-1), Cauchy-Binet gives

det H_q = sum_{S, |S|=q} [prod_(r in S) a_r]
          [prod_(r<s in S) (lambda_r-lambda_s)^2].

For exactly q visible poles this reduces to product(residues) times the squared
Vandermonde.  Thus inverse conditioning necessarily degenerates under either:
(1) weak projected residue, or (2) near-coincident poles.  This is independent
of the separate static-tree short-edge margin.

## Two-pole Prony map

For two poles, s1=lambda1+lambda2 and s2=lambda1 lambda2 obey

[m1 -m0; m2 -m1] [s1;s2] = [m2;m3].

The determinant of this system is exactly D_2.  Hence perturbation amplification
contains 1/D_2.  Recovering the roots of t^2-s1 t+s2 adds a further local
1/|lambda1-lambda2| factor.

## Three independent stability margins

A proof-grade dynamic inverse certificate should therefore report separately:
- topology margin: shortest internal static resistance w_min versus the
  effective-resistance perturbation bound;
- spectral separation: min_{r!=s}|lambda_r-lambda_s|;
- visibility/residue margin: projected residue weights or matrix singular-value
  lower bounds.

No one of these substitutes for the others.

## Storage intervals

H_0=sum_r R_r=L_BI C_I^{-1} L_IB.  If boundary leaf i has conductance g_i to a
single internal parent v, then (H_0)_ii=g_i^2/c_v.  Therefore if
|Hhat_0(ii)-H_0(ii)|<=rho_H, g_i in [g_-,g_+], and Hhat_0(ii)>rho_H,

c_v in [ g_-^2/(Hhat_0(ii)+rho_H),
         g_+^2/(Hhat_0(ii)-rho_H) ].

Pole-location errors do not enter H_0 if residues are independently bounded;
they enter H_1=sum lambda_r R_r and all higher stripping moments.

## Zero mode

The required zero mode is the static one Lambda(0) 1=0.  With positive internal
storage, Lambda(z)1 generally differs from zero for z !=0 because a uniform
boundary drive charges internal storage.  Imposing a dynamic zero mode is an
additional model assumption.
