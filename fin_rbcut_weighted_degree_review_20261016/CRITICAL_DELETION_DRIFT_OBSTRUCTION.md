# New bounded obstruction: shell deletion is not critically drift-close

Date: 2026-10-16. Analytical stationary model corollary; not a new FIN law,
limiting process theorem or IBA counterfamily.

## 1. A correctly supported comparison

Take the native reflected rates on the permitted opposite-parity H edges,
delete forbidden F edges, and set **all same-parity rates to zero**. Call
this L†. This avoids the support mistake of retaining phase-crossing
same-parity edges. It is a well-defined finite reversible comparison with
the same pi_D and total rate <=N; it need not be irreducible.
Let d†=L†m, d0=L0_D m, and G†=N Gamma†(m,m).

The shell theorem gives

    E_pi|d†-d0|<=2 pi(Shell_N)->0,
    0<=G0-G†<=4 1_Shell_N.

Native G0->2 in stationary L1 is inherited from UF/KR, so G†->2 in L1 too.
Neither fact gives exact native parity drift or a critical drift theorem.
The Fourier shell is not assumed exponentially rare.

## 2. Phase stationarity defects have a different, exponential bound

Every phase-changing transition of L† has an endpoint in the actual
departure layer E_N. Reversibility and the total-rate cap give stationary
phase-crossing intensity <=2N pi(E_N). Since |m|<=1,

    |N E[1_A L†m^2]|=|N E[m^2 L†1_A]|
       <=2N^2 pi(E_N)->0.

The actual MT departure estimate pays this exponential term. In particular,
using `L†m²=2m d†+Gamma†(m,m)` and the L1 bracket limit gives

    N E[1_A m d†]->-a,   N E[1_B m d†]->-b,           (1)

where a=1/2+mu/8, b=1-a, mu=E_rho0 exp(-Z²).
No conditional invariance of a leaking phase is assumed.

## 3. Compare with the already paid native drift

MT moment UI, the native weighted drift remainder and the fixed phase limit
pay `N m d0=Z_N b0(Z_N)+o_L1(1)`, b0(z)=-z³/3. Define
`M=E_rho0[Z² exp(-Z²)]>0`, rho0 proportional to exp(-z⁴/12).
Integrating `(z exp(-z²)rho0)'` gives

    E[(Z⁴/3-1)exp(-Z²)]=-2M.

With phase probability p(z)=1/2+exp(-z²)/8, it follows that

    -2N E[1_A m d0]->2a-M/2,
    -2N E[1_B m d0]->2b+M/2.                          (2)

Subtracting (2) from (1) proves the signed limits

    N E[1_A m(d†-d0)]->-M/4,
    N E[1_B m(d†-d0)]->+M/4.

Consequently

    liminf E_pi | Z_N N^(3/4)(d†-d0) | >= M/2 >0.    (3)

These are stationary weighted-drift statements on the actual local count
family. They are not inferred by multiplying a TV error by N. All expectations
here are under the exact stationary local law; comparison to iid-start
parent histories is not claimed for L†.

## 4. What (3) rules out, and what it does not

Unscaled L1 drift closeness and vanishing deleted shell mass do not justify
critical drift transfer. The simple deletion construction fails the exact
KR drift gate and even a natural **weighted critical L1** approximation gate,
by a fixed positive amount. No shell anti-concentration rate can remove the
signed static mismatch while preserving all the preceding premises.

An admissible exact repair must therefore compensate the drift, not simply
discard a rare Fourier interface. It must also change the bracket on bulk
mass if it realizes the intended near-closed phases, as the older BA/FP
conditional calculation already shows. This does not prove that every
possible repair is impossible, produce different limiting SDEs, refute
native HB, establish IBA failure or derive a physical memory device.
