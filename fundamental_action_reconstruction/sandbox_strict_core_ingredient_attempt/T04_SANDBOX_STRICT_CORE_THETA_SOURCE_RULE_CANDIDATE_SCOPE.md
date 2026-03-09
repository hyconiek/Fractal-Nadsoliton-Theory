# T04 Sandbox Strict-Core Theta-Source Rule Candidate Scope

Status: `T04_SANDBOX_STRICT_CORE_THETA_SOURCE_RULE_CANDIDATE_SCOPE_NO_FALSE_PASS`
As of: `2026-03-09`

## Goal

After `F03/P03/N03`, the sandbox already contains one non-placeholder
strict-core theta-source skeleton attempt.

The next narrower question is:

```text
can the sandbox write one conditional strict-core theta-source rule candidate
without importing the axiom-augmented branch and without contradicting C50?
```

This step does **not** try to export actual `theta_1`, `theta_2`.
It only tries to sharpen the strict-core-only skeleton into one conditional
rule form.

## Support reused

1. `C33`
   - local phase formula class
     `theta_i = atan2(<s_i,u_i>, <c_i,u_i>)`,
2. `C34`
   - local representative class
     `u_i(theta_i) = cos(theta_i)c_i + sin(theta_i)s_i`,
3. `C49`
   - packet-ready conditional populated-instance schema,
4. `C50`
   - denial of a packet-ready strict-core minimal source skeleton for actual
     `theta_1`, `theta_2`,
5. `F03`
   - current sandbox non-placeholder skeleton attempt.

## Question

Is the following narrow move honest?

```text
if a strict-core populated basis-pair instance (u_1,u_2) were supplied,
then theta_1, theta_2 would be serialized by the C33 atan2 formulas
```

This is weaker than a source discharge.
It is only one conditional rule candidate above the strict-core skeleton
attempt.

## Hard limits

`T04` must not claim:

1. actual `theta_1`, `theta_2`,
2. actual populated `u_1`, `u_2`,
3. actual internal orientation datum,
4. actual `E_orient`,
5. admissible `S_sel_int`,
6. strict-core selector closure,
7. ToE closure.

It must also avoid importing the axiom-augmented phase branch.
