#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2669_s1619_boundary_cycle_cross_invariant_boolean_ansatz_audit.json"
MD = GEN / "p2669_s1619_boundary_cycle_cross_invariant_boolean_ansatz_audit.md"
P2668 = GEN / "p2668_s1618_orientation_datum_to_boundary_cycle_cross_invariant_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NEGATIVE_EXPORT_FLAGS = [
    "boundary_cycle_boolean_cross_invariant_exported",
    "orientation_datum_to_boundary_cycle_cross_invariant_exported",
    "canonical_pair12_boundary_orientation_map_exported",
    "orientation_reversal_forbidden_internally",
    "boundary_phase_bit_target_exported_unconditionally",
    "intrinsic_entropy_level_exported",
    "bit_to_action_map_sourced_unconditionally",
    "bit_to_length_map_sourced_unconditionally",
    "target_independent_beta_source_exported",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "toe_closure_claimed",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, ".",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json",
        "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
    ], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:120]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "cross_invariant_boolean_content": "cross-invariant|Boolean|GF\\(2\\)|ANF|mixed invariant|boolean ansatz",
        "boundary_cycle_sector_content": "boundary-cycle|boundary cycle|square holonomy|sector swap|boundary sector|non-exact sector",
        "pair2_selector_content": "pair2.*positive|pair1/pair2|pair12|w_break|selector witness",
        "source_origin_guard_content": "physical origin|not a convention|declared|source theorem|current strict theorem",
        "nonclosure_guard_content": "QW-2191|role-bearing L_total|ToE closure|beta source|nonexport|future-route",
    }
    return {"tool": "rg", "mode": "content-first boundary-cycle Boolean cross-invariant ansatz audit", "patterns": {key: rg_count(pattern) for key, pattern in patterns.items()}}


def upstream_consistency() -> dict[str, Any]:
    p2668 = load_json(P2668)
    decision = p2668.get("closure_decision", {})
    return {
        "p2668_existing_orientation_export_present": decision.get("existing_orientation_export_present") is True,
        "p2668_no_source_forbids_sector_swap": decision.get("source_forbids_sector_swap") is False,
        "p2668_no_source_ties_pair2_to_sector1": decision.get("source_ties_pair2_to_sector1") is False,
        "p2668_no_boundary_phase_bit_target": decision.get("boundary_phase_bit_target_exported_now") is False,
    }


def eval_anf(coeffs: dict[str, int], p: int, s: int) -> int:
    return (coeffs["c"] ^ (coeffs["p"] & p) ^ (coeffs["s"] & s) ^ (coeffs["ps"] & p & s)) & 1


def anf_terms(coeffs: dict[str, int]) -> list[str]:
    labels = {"c": "1", "p": "pair2", "s": "sector1", "ps": "pair2*sector1"}
    return [labels[key] for key, value in coeffs.items() if value]


def boolean_ansatz_witness() -> dict[str, Any]:
    rows = []
    for bits in itertools.product([0, 1], repeat=4):
        coeffs = dict(zip(("c", "p", "s", "ps"), bits, strict=True))
        truth = {(p, s): eval_anf(coeffs, p, s) for p in (0, 1) for s in (0, 1)}
        sector_swap_odd = all(truth[(p, s)] ^ truth[(p, 1 - s)] == 1 for p in (0, 1) for s in (0, 1))
        ties_pair2_to_sector1 = truth[(1, 1)] == 1 and truth[(1, 0)] == 0
        distinguishes_pair_branch = truth[(0, 1)] != truth[(1, 1)] or truth[(0, 0)] != truth[(1, 0)]
        contains_mixed_term = coeffs["ps"] == 1
        convention_free_source_exported = False
        passes = sector_swap_odd and ties_pair2_to_sector1 and distinguishes_pair_branch and contains_mixed_term and convention_free_source_exported
        rows.append({
            "coefficients": coeffs,
            "anf_terms": anf_terms(coeffs),
            "truth_table": {f"p{p}_s{s}": val for (p, s), val in truth.items()},
            "sector_swap_odd": sector_swap_odd,
            "ties_pair2_positive_to_sector1": ties_pair2_to_sector1,
            "distinguishes_pair_branch": distinguishes_pair_branch,
            "contains_mixed_pair_sector_term": contains_mixed_term,
            "convention_free_source_exported": convention_free_source_exported,
            "passes_cross_invariant_acceptance_now": passes,
        })
    mathematical_candidates = [row for row in rows if row["sector_swap_odd"] and row["ties_pair2_positive_to_sector1"]]
    mixed_candidates = [row for row in mathematical_candidates if row["contains_mixed_pair_sector_term"]]
    return {
        "statement": "P2669 enumerates the full GF(2) Boolean ansatz f(pair2, sector1)=c+a*p+b*s+d*p*s for a possible boundary-cycle cross-invariant.  Mathematical functions that are odd under sector swap and tie pair2 to sector 1 exist, but the strict sector-swap-odd condition actually eliminates the mixed p*s term in this degree-2 ansatz.  The branch-sensitive survivor uses pair and sector labels additively, so the missing object is a convention-free physical origin for those labels, not a sourced mixed term.",
        "rows": rows,
        "total_boolean_ansatz_count": len(rows),
        "mathematical_sector_swap_odd_tie_count": len(mathematical_candidates),
        "mixed_pair_sector_candidate_count": len(mixed_candidates),
        "passing_cross_invariant_count": sum(1 for row in rows if row["passes_cross_invariant_acceptance_now"]),
        "mathematical_candidates_exist": bool(mathematical_candidates),
        "mixed_candidates_exist": bool(mixed_candidates),
        "convention_free_source_exported_for_any_candidate": any(row["convention_free_source_exported"] for row in rows),
    }


def closure_decision(witness: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "P2669_BOUNDARY_CYCLE_BOOLEAN_CROSS_INVARIANT_ANSATZ_AUDIT__MATHEMATICAL_CANDIDATES_NO_SOURCE",
        "professorial_verdict": "P2669 constructs the missing object at the finite Boolean ansatz level rather than only saying it is absent.  The exhaustive GF(2) ANF audit finds mathematical cross-invariant candidates that are odd under sector swap and tie pair2 to sector 1.  However, the sector-swap-odd constraint excludes mixed p*s terms in this finite ansatz; the only branch-sensitive survivor is additive in pair and sector labels, whose physical coding origin is not exported.  Therefore the result is a precise theorem target, not a source theorem, and it exports no boundary-phase bit target, UV unit, beta source, QW-2191 discharge, L_total reopening, or ToE closure.",
        "next_honest_step": "Attempt a physical-origin theorem for the branch-sensitive additive Boolean candidate `1 + pair2 + sector1`, or prove that no higher-order invariant can do better.  The proof must derive the pair variable and boundary-sector variable from bridge-completed nadsoliton boundary dynamics, not from labels.  If that origin cannot be supplied, promote P2669 to a finite no-go: Boolean cross-invariants exist mathematically but cannot source entropy-bit anchoring without an imported coding convention.",
        "mathematical_candidates_exist": witness["mathematical_candidates_exist"],
        "mixed_candidates_exist": witness["mixed_candidates_exist"],
        "passing_cross_invariant_count": witness["passing_cross_invariant_count"],
        "convention_free_source_exported_for_any_candidate": witness["convention_free_source_exported_for_any_candidate"],
        "boundary_phase_bit_target_exported_now": False,
        "unconditional_uv_unit_selected_now": False,
        "beta_source_exported_now": False,
        "qw2191_discharged_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    witness = payload["boolean_ansatz_witness"]
    decision = payload["closure_decision"]
    lines = [
        "# P2669/S1619 boundary-cycle Boolean cross-invariant ansatz audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first repo grep",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines += ["", "## Upstream consistency"]
    for key, value in payload["upstream_consistency"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines += [
        "",
        "## Boolean ansatz witness",
        witness["statement"],
        f"Total Boolean ansatz count: `{witness['total_boolean_ansatz_count']}`.",
        f"Sector-swap-odd tie candidates: `{witness['mathematical_sector_swap_odd_tie_count']}`.",
        f"Mixed pair-sector candidates: `{witness['mixed_pair_sector_candidate_count']}`.",
        f"Passing cross-invariant count: `{witness['passing_cross_invariant_count']}`.",
        f"Convention-free source exported for any candidate? `{witness['convention_free_source_exported_for_any_candidate']}`.",
        "",
        "## Candidate rows satisfying mathematical sector-swap/tie conditions",
        "| ANF terms | truth table | mixed term? | source exported? | passes? |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    for row in witness["rows"]:
        if row["sector_swap_odd"] and row["ties_pair2_positive_to_sector1"]:
            lines.append(f"| `{row['anf_terms']}` | `{row['truth_table']}` | `{row['contains_mixed_pair_sector_term']}` | `{row['convention_free_source_exported']}` | `{row['passes_cross_invariant_acceptance_now']}` |")
    lines += [
        "",
        "## Verdict",
        decision["professorial_verdict"],
        f"Decision: `{decision['decision']}`.",
        f"Boundary-phase bit target exported now? `{decision['boundary_phase_bit_target_exported_now']}`.",
        f"Beta source exported now? `{decision['beta_source_exported_now']}`.",
        f"QW-2191 discharged now? `{decision['qw2191_discharged_now']}`.",
        f"Role-bearing L_total now? `{decision['role_bearing_ltotal_now']}`.",
        f"ToE closure now? `{decision['toe_closure_now']}`.",
        "",
        "## Next honest step",
        decision["next_honest_step"],
        "",
        "## Negative exports",
    ]
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    audit = semantic_rg_audit()
    upstream = upstream_consistency()
    witness = boolean_ansatz_witness()
    decision = closure_decision(witness)
    payload: dict[str, Any] = {
        "status": "P2669_BOUNDARY_CYCLE_BOOLEAN_CROSS_INVARIANT_ANSATZ_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": audit,
        "source_hashes": {"P2668": sha256_file(P2668), "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET), "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT)},
        "upstream_consistency": upstream,
        "boolean_ansatz_witness": witness,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2669/S1619 boundary-cycle Boolean cross-invariant ansatz guard", "## P2669/S1619 boundary-cycle Boolean cross-invariant ansatz guard\n\n`P2669/S1619` constructs the missing P2668 object at finite Boolean-ansatz level: over `GF(2)`, mathematical functions `f(pair2, sector1)=c+a*p+b*s+d*p*s` exist that are odd under sector swap and tie `pair2` to boundary sector `1`.  This is not yet a source theorem: the strict oddness constraint leaves only sector-label and additive pair/sector-label candidates, and their coding still lacks a bridge-completed physical origin.  Therefore no boundary-phase bit target, intrinsic UV unit, beta source, `QW-2191` discharge, bridge completion, role transfer, role-bearing `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2669/S1619 Boolean cross-invariant Ltotal guard", "## P2669/S1619 Boolean cross-invariant Ltotal guard\n\n`P2669/S1619` keeps `L_total` closed to Boolean cross-invariant entropy terms.  A future variational coefficient must first derive the pair variable and boundary-sector variable from nadsoliton boundary dynamics rather than from labels or coding convention.\n")
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
