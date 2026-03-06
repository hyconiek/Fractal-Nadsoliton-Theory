# O1 Explicit A1 Operator Definition Spec

## Goal

Specify the minimal missing operator-level object that would make coefficient
extraction for the `H` lane computable.

## Motivation

`H28` establishes that the repository already has:

- provenance for `Route A`
- coefficient semantics for `A_1`
- upstream targets such as `A1_cc = (A_1)_{c_1 c_1}`

but still lacks the actual operator definition from which these quantities could be
 evaluated.

## Minimal Required Object

Define an explicit operator object for `pair1`:

- `A_1_ext`

with the following contract:

1. domain/codomain:
   - acts on `V_1 = span{c_1, s_1}`
2. lane discipline:
   - `hypothesis_extension_only`
3. strict-core discipline:
   - `base_kernel_contains_obs = false`
   - `strict_core_reinterpretation = false`
4. selector discipline:
   - `selector_smuggling = false`

## Two Admissible Construction Modes

`A_1_ext` may be given only in one of two admissible ways:

### Mode A

Direct composite export:

- `A_1_ext = exported_composite_A_1`

with an explicit, finite, persisted definition.

### Mode B

Factored operator chain:

- `A_1_ext = P_1 E_1^* G_light^* R_mat^* O_readout R_mat G_light E_1 P_1`

where:

- `O_readout` is the preferred name for the internal readout/backreaction operator
- legacy alias:
  - `O_obs`

## Computability Requirement

Whichever admissible mode is chosen, the definition must support evaluation of:

- `a_1 := <c_1, A_1_ext c_1>`
- `b_1 := <c_1, A_1_ext s_1>`
- `d_1 := <s_1, A_1_ext s_1>`

and therefore:

- `trace_A_1 = a_1 + d_1`
- `Delta_1 = (a_1 - d_1, b_1)`

## Minimum Acceptance Conditions

The explicit operator definition is acceptable only if:

1. it is persisted in the repository,
2. it is finite and inspectable,
3. it does not inject `theta_1 = theta_2 = 0` by definition,
4. it allows coefficient evaluation on `pair1`,
5. it preserves the `lambda_obs -> 0` decoupling limit when the factored lane is
   used.

## Result

Current missing object is no longer vague.

It is precisely:

- an explicit, persisted, coefficient-computable operator definition for `A_1_ext`
  on `V_1 = span{c_1, s_1}`.

## Honest Frontier

- `O1_B1 := the minimal explicit operator definition needed to compute a_1, b_1, d_1 is now specified, but no persisted computable A_1_ext instance exists yet in either admissible mode`
- `H28_B1 := the current repository state contains no computable operator-level source from which a_1, b_1, d_1 can be actually exported or evaluated for pair1, even though Route A provenance and coefficient semantics are already in place`
- `T12_B1 := the typing judgment with totality and uniqueness is specified but not discharged for the current selector track`
- `T2_B1 := the bridge theorem is specified but not discharged; strict-core target slot and equivalence/export map remain absent`
- `C32_B2 := raw cross-pair overlap route is formally degenerate under the strict orthonormal-disjoint mode scaffold and thus does not export alpha_12`

## Negative Claims Maintained

- no theorem-level PASS
- no full-closure PASS
- no claim that `A_1_ext` already exists
- no claim that `O_readout` is already derived from strict core
- no claim that current `K_obs` is validated
