# EXACT-K4xK6-FACTORIZED-FAMILY-101
## The k4 crossing on the pure-k6 trunk factorizes exactly into two scalar mean-field equations

Date: 2026-09-26

Status:
- exact invariant-family reduction;
- exact scalar stationary equations;
- high-g support limits follow analytically from the factorized field;
- transverse index changes away from the factorized subspace are numerical unless stated otherwise.

## 1. Factorized field

Take

    h_j
      = J4 cos(2*pi*j/3)
        + J6 (-1)^j.

The mod-3 coordinate and parity coordinate are independent on the twelve-label
carrier. Therefore the partition function factorizes exactly:

    Z(J4,J6)
      = Z4(J4) Z6(J6).

Consequently the two stationary equations decouple.

Define

    t = 3 J4 / 2.

Then

    t
      = c4 R(t),

    R(t)
      = [exp(t)-1]/[exp(t)+2],

    c4 = g lambda4 / 4,

and independently

    J6
      = c6 tanh(J6),

    c6 = g lambda6 / 12.

So the full k4+k6 stationary family is the Cartesian product of the pure-k4
and pure-k6 scalar solution sets.

## 2. Exact k4 crossing on the k6 trunk

The small k4 root collides with t=0 at

    boxed:
    g4=12/lambda4
      =5.455614632675739.

At the same g the k6 coordinate is already nonzero.

The critical k4 plane is two-dimensional and carries the D3-type representation
inside the D6 parent.

The cubic coefficient along a normalized k4 direction is nonzero:

    T3≈-0.0554906067896.

Hence the local k4 daughters have linear amplitude, not square-root amplitude.

## 3. Two global k4 sheets for g>g4

For g>g4 the scalar k4 equation has:
- one small negative root;
- one large positive root;
- the zero root.

After multiplying by the nonzero k6 solution, these give three distinct
stationary families.

### Negative-k4 sheet

At large g it concentrates on

    S4={2,4,8,10}

with

    p*=1/4 on each label.

Therefore

    asymptotic Morse index = 3.

Its D12 stabilizer has order 4, so the orbit size is 6.

### Positive-large-k4 sheet

At large g it concentrates on

    S2={0,6}

with

    p*=(1/2,1/2).

Therefore

    asymptotic Morse index = 1.

Its stabilizer also has order 4 and orbit size 6.

### Zero-k4 sheet

This is simply the pure-k6 trunk, with large-g support

    S6={0,2,4,6,8,10}

and asymptotic index 5.

## 4. Why one local crossing feeds very different support complexities

The D3 cubic does not uniquely determine the global support destination.

The two nonzero scalar k4 sheets have opposite signs.

For positive J4, the field selects one mod-3 class; after intersecting with the
positive-parity k6 support, only two labels survive.

For negative J4, two mod-3 classes are selected; after the same parity
intersection, four labels survive.

Thus the exact factorization explains algebraically why the same local
representation can feed either

    support size 2

or

    support size 4.

## 5. Relation to previous branches

The positive-large sheet is exactly the large-g d=6 pair branch later used in
report 92.

The negative sheet is a previously unresolved four-label parent. Its own
transverse breaking is classified in report 102.

So the k4 crossing is a genuine branch splitter inside the high-g k6 trunk.
