# N501 Current First Strict `R1` Target‑Slot Span Sign Gauge‑Irrelevance Theorem

Status: `N501_DISCHARGED_CURRENT_FIRST_STRICT_R1_TARGET_SLOT_SPAN_SIGN_GAUGE_IRRELEVANCE_THEOREM_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

After `B2/B3`, the remaining strict uniqueness gap is narrower:

- continuous `O(2)` axis choice is already cut to residual `Z2` (sign) on two independent strict lanes (`N492`, `N500`),
- but **sign‑sensitive physical orientation** (lifting residual `Z2`) remains open.

This theorem closes one potential ambiguity at the strict target‑slot level:

```text
for the declared residual orientation datum target slot (R1),
residual Z2 sign flips of the representative vectors do not change the target-slot object,
because the target-slot semantics is span/projector-based.
```

So, in the declared `R1` scope, residual sign is a tracked convention layer and is **not** an additional missing ingredient
for populating the target slot as a subspace object.

This theorem does **not** claim strict-core selector closure, does **not** discharge `QW-2191` globally, and does **not** claim ToE closure.

## Strict-admissible inputs reused

1. `R1`
   - strict-core residual orientation datum target-slot export packet:
     the target object class is explicitly `S_orient_cand(theta_1,theta_2)=span{u_1(theta_1),u_2(theta_2)}`,
2. `A10`
   - anti-overclaim boundary.

Probe-level consistency evidence (not needed for the proof):

- `P456` audits that the exported sigma-int theta-pair vectors are orthonormal and that the resulting span projector is
  numerically invariant under sign flips on the current exported instance.

## Setup (R1 target-slot semantics)

By `R1`, for angles `(theta_1, theta_2)` the target-slot object is the 2D subspace:

$$
S(\theta_1,\theta_2) := \operatorname{span}\{u_1(\theta_1),u_2(\theta_2)\}.
$$

The residual sign freedom is the independent choice of signs:

$$
u_1 \mapsto s_1 u_1,\quad u_2 \mapsto s_2 u_2,\qquad s_1,s_2\in\{+1,-1\}.
$$

## Theorem (sign gauge-irrelevance for the R1 span target slot)

For any nonzero vectors $u_1,u_2$ and any $s_1,s_2\in\{+1,-1\}$:

$$
\operatorname{span}\{s_1 u_1, s_2 u_2\} = \operatorname{span}\{u_1,u_2\}.
$$

**Proof.** Multiplication by a nonzero scalar does not change the span of a vector.
Therefore each generator set spans the same subspace. ∎

Equivalently, if one encodes the target slot by its orthogonal projector $P_S$, then $P_S$ is unchanged under these sign flips.

## Meaning (no false pass)

This theorem means only:

1. in the declared `R1` target-slot semantics, residual `Z2` sign flips of the representative vectors do not change the target-slot object,
2. therefore, lifting residual `Z2` to a **sign-sensitive physical orientation** is not required for the narrow task “populate the `R1` target slot as a span object”,
3. any downstream claim that truly requires a sign-sensitive directed orientation must still either:
   - prove sign gauge-irrelevance for that downstream observable, or
   - derive an additional strict sign-sensitive datum (cf. `B3`),
   before promotion.

## What N501 does not claim

`N501` does not claim:

1. sign-sensitive physical orientation datum derived (it explicitly remains open),
2. strict-core selector closure / admissible `S_sel_int`,
3. global discharge of `QW-2191`,
4. ToE closure.

