# P1971 Strict Higgs Branch-Lock / Invertibility Lemma Packet

Status: `CONDITIONAL_HIGGS_BRANCH_LOCK_LEMMA_PROVED__POSITIVITY_PREMISE_NOT_EXPORTED__PO2_STILL_OPEN`  
As of: `2026-05-18`

## Goal

Continue from `P1970`, which identified the remaining `PO2` leftover as a
Higgs-background branch difference `delta_H` and a Yukawa-background term
`y_f*delta_H`.

`P1971` asks the next honest question:

```text
Can the current strict/PO3 exports force delta_H = 0?
```

## Conditional branch-lock lemma

For the frozen Higgs force map

```text
F(H;R) = mu_H_sq*H + lambda_H*H^3 + xi_H*R*H,
```

the script proves the identity

```text
F(H_A;R) - F(H_B;R)
= (H_A-H_B) * ∫_0^1 F'(H_B+t(H_A-H_B);R) dt.
```

Therefore, if the admissible branch interval exports a positive lower bound

```text
F'(H_B+t(H_A-H_B);R) >= m_eff_min > 0
```

for all `t in [0,1]`, then equal Higgs EOM rows imply:

```text
H_A - H_B = 0.
```

This would remove the `P1970` leftover `y_f*delta_H`.

## Current-export audit

`P1971` checks the current PO3 exports.  The result is:

```text
PO3 formal domain nonempty = true
Higgs positive-derivative / invertibility premise exported = false
```

So the branch-lock lemma is mathematically available, but its required premise
is not yet an exported strict-domain fact.

## ToE progress

This is real progress because the vague `Omega_unexported` from `P1965` is now
sharpened into a precise missing object:

```text
positive Higgs-row derivative / branch-lock invertibility on the admissible
PO3 branch interval.
```

But it is not full `PO2` closure.

## No false pass

`P1971` does not claim:

1. full `PO2` sufficiency,
2. global background-independence closure,
3. strict-core ToE closure,
4. `QW-2191` discharge.

## Output artifacts

- Script:
  `p1971_s921_strict_higgs_branch_lock_invertibility_lemma_probe.py`
- Witness:
  `generated/p1971_s921_strict_higgs_branch_lock_invertibility_lemma_probe.json`

## Next honest step

Build `P1972`: attempt to export a concrete `m_eff_min > 0` witness from the
strict coefficient/domain data.  If that fails, freeze `delta_H` as the exact
`PO2` normal-form obstruction term.

## Lay explanation

We proved what kind of lock would force the two Higgs backgrounds to be the
same: the Higgs equation must be one-to-one on the allowed range.  But the repo
has not yet proved that this one-to-one condition actually holds on its allowed
strict domain.  So the missing piece is now very precise, but still missing.
