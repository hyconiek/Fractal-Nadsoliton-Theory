# CONSERVATIVE-HYDRODYNAMIC-FDT-261
## The closed SWAP network yields two conserved density fields with diffusion and a current-noise matrix fixed by the same categorical covariance

Date: 2026-09-27

Status:
exact microscopic current calculation and continuum hydrodynamic identification.

Let

    p_a(i,t)

be the local probability / coarse density of color a in Z3.

Pure SWAP on each edge at rate rho/2 gives EXACTLY at the one-point level:

    d p_a(i)/dt
      =
      D[
        p_a(i+1)+p_a(i-1)-2p_a(i)
      ],

with

    boxed:
    D=rho/2.

In the continuum lattice-scale limit:

    boxed:
    partial_t p_a
      =
      D partial_x^2 p_a.

Because

    p_0+p_1+p_2=1,

there are exactly TWO independent conserved density fields.

## 1. Microscopic current covariance

For one edge event with left color x and right color y, define the color-current vector from left to right:

    J
      =
      e_x-e_y.

Under local equilibrium with color distribution p:

    E[J J^T]
      =
      2[
        diag(p)-p p^T
      ].

Multiplying by the edge event rate rho/2 gives the current covariance per unit time:

    boxed:
    rho[
      diag(p)-p p^T
    ]

      =
    2D S(p),

where

    S(p)
      =
      diag(p)-p p^T.

Direct enumeration at:
- p=(1/3,1/3,1/3);
- p=(0.6,0.3,0.1);
- p=(0.8,0.1,0.1);

agrees to Frobenius residual below 1.2e-16.

## 2. Fluctuating hydrodynamics

The natural conservative stochastic field equation is therefore

    partial_t p_a
      =
      -partial_x j_a,

    j_a
      =
      -D partial_x p_a
      +
      xi_a,

with

    boxed:
    <xi_a(x,t) xi_b(x',t')>
      =
      2D S_ab(p)
      delta(x-x')
      delta(t-t').

At uniform p=(1/3,1/3,1/3):

    S
      =
      diag(1/3)-11^T/9,

with eigenvalues

    0,
    1/3,
    1/3.

The zero direction is the fixed total density.
The two tangent directions are the two hydrodynamic fields.

## 3. Structural FIN connection

The same categorical covariance

    S(p)=diag(p)-p p^T

already appeared in the one-unit Onsager/FDT construction.

Here it reappears as the conservative current mobility/noise matrix.

This is a genuine structural continuity between:
- open local relaxation;
- closed conservative transport.

It is still a classical dissipative field theory, not quantum field theory.
