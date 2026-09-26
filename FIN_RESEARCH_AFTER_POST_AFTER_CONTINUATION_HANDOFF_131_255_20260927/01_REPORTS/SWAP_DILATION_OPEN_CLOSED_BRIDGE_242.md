# SWAP-DILATION-OPEN-CLOSED-BRIDGE-242
## The same local reversible gate interpolates between open MaxEnt refresh and closed conservative transport

Date: 2026-09-26

Report 241 is not an unrelated new kinetic law.

It is the closed-environment version of the existing open Z3 refresh.

## 1. Fresh environment

Take a system trit x and a fresh environment trit e drawn uniformly.

Apply the unique swap gate:

    (x,e)
      ->
    (e,x).

After discarding the environment output, the new system state is exactly
uniform and independent of x.

Therefore a Poisson stream of fresh uniform ancillas at rate rho gives EXACTLY

    Q3=rho(U-I).

The reset-kernel replay has zero numerical residual.

## 2. Reuse the environment instead of discarding it

If the ancilla trits are retained as other units and swaps occur locally, no
information is erased.

Instead:
- old labels move into neighboring degrees of freedom;
- local relaxation becomes transport;
- finite propagation geometry matters;
- long wavelength modes appear.

Thus the architecture is

    open subsystem:
      swap + fresh ancilla
        -> MaxEnt heat bath

    closed multiunit system:
      swap + retained neighboring ancilla
        -> diffusion / memory / conserved information.

This directly implements the working FIN intuition:

    information apparently lost by a subsystem
      is transferred into relations/environment.

## 3. New conserved quantities

Pure swaps preserve the exact global counts

    N0,N1,N2,

with

    N0+N1+N2=n.

Hence the closed n-site system has

    binomial(n+2,2)

count sectors.

This conservation was NOT visible in the autonomous open Q3 generator.

It appears only when the environment is retained explicitly.

## 4. Interpretation boundary

The count conservation is both:
- the mechanism that allows long-wavelength transport;
- a new global structural prediction of this particular reversible closure.

FIN has not yet proved that this minimal swap closure, rather than a larger
reversible environment gate, is physically fundamental.

So the result is a strong candidate bridge, not final ontology.
