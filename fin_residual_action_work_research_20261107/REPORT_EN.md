# Residual action work: a conditional critical-window damping law

Date: 2026-11-07. Parent: 400c26f5d7cfe4703797d900d636ba35b5e9ff14.

The remaining work terms D_I,D_J in RW are bounded by
C_S epsilon^3[1+log(1/epsilon)] for 0<=t<=S epsilon^-2.
This is o(epsilon^2), the scale needed to transfer the three-channel
radiation polynomial to the actual corrected actions of the supplied RP
model. This is a self-audited conditional analytical result; independent
proof review remains outstanding.

## What changed

The previous absolute estimate was only O(epsilon^2). The new proof moves
the full linear propagator onto the localized work profiles. Their
modulated primitives have a uniform spatial Linfinity bound growing at
most logarithmically. Actual modal envelopes have amplitude and total
variation O(epsilon^2), so the backward work tests cost only
epsilon^2 log(2+T).

The complete residual source, after the exact continuum-input quadratic
normal form, has time-integrated L1 source/preimage norm O(epsilon).
This includes chart reconstruction, mixed bound/continuum terms,
physical cubic/quintic forces, localized commutators and both odd
normal-form channels. Physical projections are retained; free mass
multipliers are not commuted with them. The continuum correction starts
at zero and contributes O(epsilon^3) directly.

The actual actions consequently approximate

    i'=-k_a i^2-k_x i j,
    j'=-k_b j^2-k_x i j

with uniform error O_S(epsilon^3[1+log(1/epsilon)]) on the critical window.
The coefficients are exactly CD's three-channel coefficients, without
refitting. For the existing preparation, the leading I profile is

    (omega epsilon^2/2)/(1+(k_a omega/2)epsilon^2 t).

No pointwise monotonicity of exact actions, waveform approximation or
all-time scattering is asserted.

## Reproducible checks and their limits

`check_work.py` passed 39 checks: 29 predecessor integrity checks, exact
residual-source identities and channel coefficients, and grouped numerical
controls. The spectral diagnostic evaluates 36 profile/carrier/time cases
at w=1.25,1.5,1.75 and times 8,32,128. It uses the exact Darboux transform
and time-integrated Fourier propagator on a finite periodic grid.

- Maximum primitive derivative discrepancy: 2.82e-9.
- Maximum sampled supremum change when doubling grid resolution: 4.89e-5.
- Maximum sampled supremum change when enlarging the domain: 1.17e-15.

These are finite-grid consistency controls, not certified uniform
constants, a full nonlinear PDE simulation, or a proof of the theorem.
The analytical source ledger and primitive estimate carry the argument.

## Scope and next step

The claim is conditional on the predecessor global estimates, exact
canonical/continuum maps and critical local radiation bound. All earlier
files and AGENTS hashes are preserved. Historical OPEN dispositions in RW
remain historical; this package records their new conditional resolution.

Next perform a focused proof audit of AW and its dependencies, especially
the two odd-channel normal forms, weighted source estimates and the
uniform-in-space primitive bound. The physical sourcing, units and
empirical FIN interpretation remain separate open obligations.
