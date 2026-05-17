# P1963 Strict PO3 Double-Run Machine Checker Packet

Status: `PO3_FORMAL_DOMAIN_NONEMPTY_MACHINE_CHECKED__GLOBAL_STRICT_CORE_CLOSURE_STILL_OPEN`
As of: `2026-05-17`

## Purpose

`P1945` required a real reproducible checker run for:

```text
THM_A_EPS_NONEMPTY_V1
```

`P1963` supplies that run.  It performs two deterministic machine checks of
the same formal admissible-branch domain:

```text
D_adm = branches satisfying strict tuple, invariant-triplet closure,
        shared scheme, and rho_i <= eps_star

A_eps = bounded admissible subclass
```

and verifies the existential witness:

```text
exists b in D_adm such that b in A_eps
```

with witness:

```text
b := BR_C_strict_consistent_seed_machine_checked_v1
```

## Machine Result

The checker exports:

```text
checker_exit_code = 0
proof_artifact_hash_sha256 = stable across run1/run2
domain_encoding_digest = stable across run1/run2
double_run_gate_pass = true
```

Therefore the PO3 restamp is:

```text
PASS_MACHINE_CHECKED_FORMAL_DOMAIN_NONEMPTY
```

## Exact Scope

This certifies only the formal encoded non-emptiness target from
`P1939-P1945`.

It does not claim:

```text
global background-independence closure
full renormalization closure
full Cutkosky unitarity
QW-2191 selector discharge
ToE closure
```

## Open Blocks After P1963

The minimum remaining strict-core open blocks are reduced only in the formal
PO3 sense:

```text
R1 renormalization theorem-grade closure
U1 unitarity/Cutkosky theorem-grade closure
B1 background-independence global theorem closure
S1 selector obstruction QW-2191 discharge
PO2 sufficiency derivation for branch-policy closure
cross-scheme finite-part transport theorem linking R1/U1/B1
```

## Outputs

- `p1963_s913_strict_po3_double_run_machine_checker.py`
- `generated/p1963_s913_strict_po3_double_run_machine_checker.json`
- `generated/p1963_s913_strict_po3_run1_meta.json`
- `generated/p1963_s913_strict_po3_run2_meta.json`
- `generated/p1963_s913_strict_po3_repro_compare.json`

## Next Honest Step

Build `P1964`: attack `PO2` sufficiency by deriving `DELTA_BG_Yf=0` from
`C1-C4` inside the full non-skeleton EOM trace, using the now machine-checked
nonempty `A_eps` domain as a witness class rather than as global closure.

## Lay Explanation

One formal blocker is now genuinely checked by a reproducible machine run:
the allowed branch class is not empty in the encoded domain.  This is not the
whole theory.  It only says that one required mathematical room has at least
one valid point inside it; the larger proofs for renormalization, unitarity,
background-independence, and selector closure remain open.
