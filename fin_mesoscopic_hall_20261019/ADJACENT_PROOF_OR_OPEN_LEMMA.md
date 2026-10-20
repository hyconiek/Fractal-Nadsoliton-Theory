# ADJACENT is paid after exact carry/source cancellation

Date: 2026-10-20. Native HB/source, not the HYP source. All normalization
factors and the reflected source s_x=(N/2)W_x L0^D m(x) are retained.

## 1. Compare conditional masses, not separate unnormalized limits

Let beta_l=W(A_ball)/W_(A,l) and beta_(l+1)=W(B_right)/W_(A,l+1).
By MH-W(4),(5), uniformly,

    beta_l=b_N(m_l)/A_N(m_l,z_l)+O(h^(4/3)),
    beta_(l+1)=b_N(m_(l+1))/A_N(m_(l+1),z_(l+1))+O(h^(4/3)).

The common normalization W_N(m,0) cancels EXACTLY within each ratio.
The smooth integral ratio has bounded m,z derivatives. Its adjacent
increments are m_(l+1)-m_l=2/N=2h^2 and z_(l+1)-z_l=2N^-3/4=2h^(3/2).
Lattice residue changes are already paid uniformly by Poisson, so

    |beta_l-beta_(l+1)|<=C h^(4/3).                 (1)

In particular this is o(h). An unrelated pair of o(1) Riemann limits
would not have proved (1).

## 2. Actual native drift source on these fibres

The UF scalar cancellation on the parity line, finite LOO error and
absence of a linear fast term give on the growing central cutoff

    |(N/2)L0 m|<=C[N|m|^3+||Y||^2+1]
                    <=C[N^(1/4)+||Y||^2+1].        (2)

The additional UF term ||w|| |m|^3 is absorbed in the same bound.
Conditional fast moments on a compact slice are uniform by the native
Gaussian envelope. Outside the cutoff use the actual cap |(N/2)L0^D m|<=N
and charge the superpolynomial weight tail; at E the exponential estimate
is retained. Thus, for the ENTIRE phase fibre,

    sum_(X_(A,l)) |s_x|<=C N^(1/4) W_(A,l).         (3)

The inherited PC compact positive-current bound gives
J_(A,l),J_(A,l-1) comparable to N W_(A,l), uniformly. This is a property
of the fixed marginal ansatz and exact native source, not a candidate
rate cap or an assumed closed phase generator.

## 3. Cancellation before any error bound

Without E redistribution the exact marginal formulas give

    O_l(A_ball)=J_(A,l-1) beta_l+s(A_ball),
    I_(l+1)(B_right)=J_(A,l) beta_(l+1),
    J_(A,l)-J_(A,l-1)=s(X_(A,l)).

Therefore their difference is EXACTLY

    J_(A,l-1)(beta_l-beta_(l+1))
        +s(A_ball)-s(X_(A,l)) beta_(l+1).           (4)

Use (1),(3) and the positive carry bound. The first term is
O(h^(4/3)J_(A,l)); the two source terms are O(N^-3/4 J_(A,l)).
The correction to incoming/outgoing mass is bounded by 2N W(E) on each
fibre, hence exponentially smaller. Consequently

    |O'_l(A_ball)-I'_(l+1)(B_right)|
       <=C[N^-2/3+N^-3/4+poly(N)exp(-bN)]J_(A,l)
        =o(N^-1/2 J_(A,l)).                        (ADJACENT)

This retains exact statewise source and true departure redistribution.
It does not compare unsupported full path laws, infer a Markov scalar
quotient, or construct transport on every cut.
