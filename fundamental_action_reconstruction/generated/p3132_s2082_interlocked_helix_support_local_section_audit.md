# P3132/S2082 Interlocked-helical support-local-section audit

Status: `P3132_INTERLOCKED_HELIX_ZETA_OS_RELATIVE_LOCK_NO_ABSOLUTE_SECTION`

## New typed candidate

`Zeta_OS^HL` models the useful part of the screw/nut intuition.  Two helical
carriers may possess a permitted relative registry and thereby block
**relative** rotation.  The question is whether that relation also provides
the missing strict source of an absolute support origin and a polarity.

This is a finite `Z12` abstraction, not a claim that DNA or an ordinary
mechanical screw is already a nadsoliton source law.

## Finite model

Let

```text
H = (p_A, a, epsilon_A; p_B, b, epsilon_B),
delta(H) = b - a mod 12,
L(p_A, p_B, delta) = locked relative registry.
```

For a fixed allowed registry `delta_0`, the locked locus is

```text
{(a, a + delta_0) : a in Z12}.
```

It has 12 states.  A common rotation/translation acts by

```text
t . (a,b) = (a+t,b+t),
```

and preserves both the registry and the lock.  Thus the interlock removes
independent relative rotation but leaves a full twelve-element diagonal
global-phase orbit.

## No-absolute-section result

Suppose a source-free quotient-level rule selected an absolute origin
`s([H])`.  Quotient invariance would give

```text
s([a,b]) = s([a+t,b+t]),
```

whereas origin covariance would require

```text
s([a+t,b+t]) = s([a,b]) + t.
```

For `t=1` these are incompatible in `Z12`.  The raw rule `s(H)=a` is not a
counterexample: it imports the phase label instead of deriving a source.

Inversion also pairs the two handedness signs unless a nonzero strict
inversion-odd source law selects one polarity.

## Decision

The screw/nut-like relation is a real new **relative-lock** object, but it is
not an absolute-origin source.  It cannot currently supply `Zeta_OS`,
`Gamma_SO`, or a discharge of `QW-2191`.

The precise missing object is now sharper: a strict, nontranslation,
support-local **helical-lock defect** `D_HL` would have to provide both a
nonzero inversion-odd value and a support representative.  A chosen axis,
boundary clamp, base label, or external mechanical environment would be an
imported premise rather than that source.

No unit-bearing action, spacetime EOM, `L_total`, bridge completion,
role transfer, or ToE closure is exported.
