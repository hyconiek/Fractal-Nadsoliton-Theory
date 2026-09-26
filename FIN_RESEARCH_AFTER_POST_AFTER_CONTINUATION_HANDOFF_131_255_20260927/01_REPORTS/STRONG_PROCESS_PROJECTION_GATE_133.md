# STRONG-PROCESS-PROJECTION-GATE-133
## PLJ=Lc is insufficient; exact self-adjoint closure and observable memory are mutually exclusive

Date: 2026-09-26

Status:
- exact linear-algebra theorem;
- exact three-state continuous-time Markov counterexample;
- supplies the corrected gate for any future task formerly called 132.

## 1. Exact Markov counterexample to the weak gate

Use the fine positive Laplacian

        [ 1   0  -1 ]
    L = [ 0   3  -3 ]
        [-2  -2   4 ]

whose rows sum to zero and whose off-diagonal entries are non-positive.

Define the coarse embedding and projection

        [1 0]
    J = [1 0],
        [0 1]

        [1/2 1/2 0]
    P = [0    0   1],

so

    P J = I.

Then

        [ 2 -2]
    Lc =[ -4 4]

satisfies exactly

    P L J = Lc.

But

    P L^2 J - Lc^2

is

        [1 -1]
        [0  0],

which is nonzero.

Therefore the projected semigroup differs already at order t^2.

Numerically,

    ||P exp(-0.1 L) J - exp(-0.1 Lc)||_2
      ≈ 0.00517082,

    at t=0.5:
      ≈ 0.0506651,

    at t=1:
      ≈ 0.0899286.

So the weak matrix condition is not a process-level closure criterion.

## 2. Correct exact criterion

A sufficient exact criterion is the intertwining

    L_fine J
      =
      J L_coarse.

Then automatically

    exp(-t L_fine) J
      =
      J exp(-t L_coarse)

for all t>=0.

For a self-adjoint positive L_fine with isometric J and P=J*,
this condition is also forced by exact projected-semigroup equality.

Indeed, write

    Pi=J J*.

The second derivative defect is

    J* L^2 J
      -(J* L J)^2

    =
      J* L(I-Pi)L J

    =
      [(I-Pi)LJ]^*
      [(I-Pi)LJ]
      >=0.

It vanishes if and only if

    (I-Pi)LJ=0,

that is, the coarse subspace is invariant.

For self-adjoint L its orthogonal complement is invariant as well.

Hence exact self-adjoint coarse semigroup preservation gives a reducing direct
sum and no dynamical feedback from the hidden complement.

## 3. Memory and exact closure are different architectures

If the coarse-hidden coupling is nonzero, block the positive generator as

        [A  B ]
    L = [B* D ].

For

    x_dot=-A x-B y,
    y_dot=-B* x-D y,

eliminating y yields

    x_dot(t)
      =
      -A x(t)
      -B exp(-Dt)y(0)
      +integral_0^t
         B exp[-D(t-s)] B*
         x(s) ds.

Thus the same coupling that breaks exact autonomous closure produces:
- preparation dependence;
- a memory kernel;
- multitime effects.

So one cannot simultaneously demand:
1. an exactly unchanged autonomous coarse semigroup; and
2. nontrivial feedback from the discarded fiber.

## 4. Corrected gate for future multiscale FIN work

Every proposed coarse-graining must choose explicitly between:

### Exact Markov closure
Prove an intertwining/lumpability condition.

### Approximate Markov closure
Give a norm/error bound over a declared time window.

### Non-Markov effective dynamics
Retain the derived memory kernel and preparation term.

The weak condition

    P L_fine J=L_coarse

alone is no longer an acceptable acceptance test.
