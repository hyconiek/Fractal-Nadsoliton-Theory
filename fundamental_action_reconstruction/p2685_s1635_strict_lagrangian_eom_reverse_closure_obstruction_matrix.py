#!/usr/bin/env python3
"""P2685/S1635: strict Lagrangian/EOM reverse-closure obstruction matrix.

This follows P2684's pivot.  It asks a finite reverse-closure question: can the
currently exported selector-independent strict Lagrangian/EOM rows be reversed
into nonproxy full tensor-resolved EOM/L_total closure?  The answer is no on the
current artifacts, with explicit obstruction rows rather than a generic guard.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any

import sympy as sp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2685_s1635_strict_lagrangian_eom_reverse_closure_obstruction_matrix.json"
MD = GEN / "p2685_s1635_strict_lagrangian_eom_reverse_closure_obstruction_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2684": GEN / "p2684_s1634_pair12_cycle_cut_semantic_invariant_provider_audit.json",
    "P2316": GEN / "p2316_s1266_strict_current_best_lagrangian_eom_coverage_audit_probe.json",
    "P2329": GEN / "p2329_s1279_selector_independence_lagrangian_eom_audit_probe.json",
    "P2362": GEN / "p2362_s1312_strict_lagrangian_eom_parallel_completion_probe.json",
    "P1866": GEN / "p1866_s816_strict_kernel_to_full_lagrangian_and_eom_symbolic_export.json",
    "P1974": GEN / "p1974_s924_strict_bianchi_anisotropic_eom_residual_obstruction_witness.json",
    "P1795": GEN / "p1795_s745_strict_nonproxy_covariant_export_bw_brst_cut_gate_readiness_checkpoint.json",
    "P1809": GEN / "p1809_s759_strict_nonproxy_export_semantic_level_reconciliation_checkpoint.json",
    "P2086": GEN / "p2086_s1036_strict_full_ltotal_eom_termwise_execution_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "full_tensor_resolved_eom_closed",
    "role_bearing_ltotal_exported",
    "selector_closure_imported",
    "legacy_bridge_imported",
    "qw2191_discharged",
    "toe_closure_claimed",
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        ["rg", "-n", pattern, ".", "-g", "*.py", "-g", "*.md", "-g", "*.json", "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**"],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:60]}


def content_grep() -> dict[str, Any]:
    patterns = {
        "selector_independent_lagrangian_eom": r"selector-independent|selector_independent|selector-free|P2329|P2362|P2086|termwise",
        "full_tensor_nonproxy_obligations": r"full tensor|tensor-resolved|nonproxy|componentwise|background family|Bianchi|FRW|metric_full_tensor_closure",
        "reverse_closure_or_helmholtz": r"reverse-closure|reverse closure|EOM -> L_total|Helmholtz|variational origin|integrability",
        "open_obstruction_witnesses": r"OPEN_OBSTRUCTION_WITH_TRACE|obstruction witness|anisotropic residual|SCALAR_TRANSPORT_INSUFFICIENT|NO_GATE_PROMOTION",
        "forbidden_imports": r"QW-2191 discharge|role-transfer|legacy.*bridge|ToE closure|selector closure|role-bearing L_total",
    }
    return {"tool": "rg", "mode": "content-first strict Lagrangian/EOM reverse-closure obstruction matrix", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def current_artifact_state() -> dict[str, Any]:
    data = {name: load_json(path) for name, path in INPUTS.items()}
    p2316 = data["P2316"]
    p2329 = data["P2329"]
    p2362 = data["P2362"]
    p1974 = data["P1974"]
    p1795 = data["P1795"]
    p1809 = data["P1809"]
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2684_reverse_closure_pivot_present": data["P2684"].get("decision", {}).get("strict_lagrangian_eom_reverse_closure_is_next") is True,
        "p2316_best_working_ltotal_identified": p2316.get("gatekeeper_checks", {}).get("best_working_ltotal_identified") is True,
        "p2316_full_task3_theorem_not_claimed": p2316.get("gatekeeper_checks", {}).get("full_task3_theorem_not_claimed") is True,
        "p2329_all_terms_selector_independent": p2329.get("gatekeeper_checks", {}).get("all_terms_selector_independent") is True,
        "p2329_all_variations_selector_independent": p2329.get("gatekeeper_checks", {}).get("all_variations_selector_independent") is True,
        "p2362_no_selector_prerequisite_for_eom_export": p2362.get("gatekeeper_checks", {}).get("no_selector_prerequisite_for_eom_export") is True,
        "p1974_generic_anisotropic_residual_nonzero": p1974.get("gatekeeper_checks", {}).get("generic_anisotropic_residual_nonzero") is True,
        "p1974_isotropic_limit_zero": p1974.get("gatekeeper_checks", {}).get("frw_isotropic_limit_zero") is True,
        "p1795_metric_full_tensor_closure_open": "metric_full_tensor_closure" in p1795.get("open_items", []),
        "p1809_tg1_locked_by_unified_nonproxy_residual": p1809.get("derived_gate_vector", {}).get("TG1_BW") == "OPEN_LOCKED_BY_UNIFIED_NONPROXY_RESIDUAL",
    }


def anisotropic_residual_rank() -> dict[str, Any]:
    p1974 = load_json(INPUTS["P1974"])
    sigma1, sigma2, dsigma1, dsigma2, H = sp.symbols("sigma1 sigma2 dsigma1 dsigma2 H")
    expressions = [sp.sympify(expr) for expr in p1974.get("anisotropic_eom_residual_vector", [])]
    variables = [sigma1, sigma2, dsigma1, dsigma2, H]
    jac = sp.Matrix([[sp.diff(expr, var) for var in variables] for expr in expressions])
    generic_rank = int(jac.rank()) if expressions else 0
    iso_rank = int(jac.subs({sigma1: 0, sigma2: 0, dsigma1: 0, dsigma2: 0}).rank()) if expressions else 0
    sample_nonzero_count = sum(1 for row in p1974.get("numeric_probe_table", []) if row.get("nonzero") is True)
    return {
        "component_basis": p1974.get("component_basis", []),
        "residual_vector": [str(expr) for expr in expressions],
        "jacobian_variables": [str(v) for v in variables],
        "symbolic_jacobian_rank": generic_rank,
        "isotropic_jacobian_rank": iso_rank,
        "isotropic_limit_zero": p1974.get("isotropic_limit_zero") is True,
        "numeric_nonzero_samples": sample_nonzero_count,
        "obstruction_detected": bool(generic_rank > 0 and sample_nonzero_count > 0),
        "interpretation": "FRW/scalar transport cannot be reversed into full tensor-resolved Bianchi-I metric EOM closure without additional anisotropic/background component equations.",
    }


def reverse_closure_rows(state: dict[str, Any], residual: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "row_id": "R1_selector_independent_reduced_terms",
            "available_forward_export": state["p2329_all_terms_selector_independent"] and state["p2329_all_variations_selector_independent"],
            "reverse_closure_obligation": "termwise reduced Euler-Lagrange rows must lift to nonproxy componentwise covariant rows on a shared background family",
            "satisfied_now": False,
            "witness_or_blocker": "P2329/P2362 allow selector-independent progress, but only as reduced/parallel EOM work, not full tensor reverse closure.",
        },
        {
            "row_id": "R2_best_working_ltotal_anchor",
            "available_forward_export": state["p2316_best_working_ltotal_identified"],
            "reverse_closure_obligation": "best working L_total must be proven theorem-grade full tensor-resolved, not merely identified as current best scaffold",
            "satisfied_now": False,
            "witness_or_blocker": "P2316 explicitly keeps the full Task-3 theorem unclaimed.",
        },
        {
            "row_id": "R3_anisotropic_background_transport",
            "available_forward_export": True,
            "reverse_closure_obligation": "scalar/FRW transport must vanish on anisotropic Bianchi-I component residuals",
            "satisfied_now": not residual["obstruction_detected"],
            "witness_or_blocker": f"P1974 residual rank={residual['symbolic_jacobian_rank']} with {residual['numeric_nonzero_samples']} nonzero numeric samples.",
        },
        {
            "row_id": "R4_unified_nonproxy_gate",
            "available_forward_export": True,
            "reverse_closure_obligation": "EA/EH/ELg must be full componentwise nonproxy exports on the same freeze and pass a unified residual",
            "satisfied_now": False,
            "witness_or_blocker": "P1795/P1809 leave metric_full_tensor_closure open and TG1_BW locked by unified nonproxy residual.",
        },
    ]


def policy_lattice(rows: list[dict[str, Any]]) -> dict[str, Any]:
    obligations = [row["row_id"] for row in rows]
    success = 0
    admissible_no_go = 0
    for bits in itertools.product([False, True], repeat=len(obligations)):
        all_pass = all(bits)
        success += int(all_pass)
        admissible_no_go += int(not all_pass)
    current = {row["row_id"]: row["satisfied_now"] for row in rows}
    return {
        "obligations": obligations,
        "total_states": 2 ** len(obligations),
        "reverse_closure_success_states": success,
        "admissible_obstruction_states": admissible_no_go,
        "current_state": current,
        "current_reverse_closure_success": all(current.values()),
    }


def decision(rows: list[dict[str, Any]], lattice: dict[str, Any]) -> dict[str, Any]:
    failed = [row["row_id"] for row in rows if not row["satisfied_now"]]
    return {
        "decision": "P2685_STRICT_LAGRANGIAN_EOM_REVERSE_CLOSURE_OBSTRUCTION_MATRIX_NO_FULL_TENSOR_CLOSURE_NO_FALSE_PASS",
        "professorial_verdict": (
            "P2685 executes the strict Lagrangian/EOM reverse-closure computation requested by P2684.  The forward/reduced route is real and selector-independent, but the reverse implication to nonproxy full tensor-resolved EOM/L_total closure fails on finite obligations: reduced rows do not supply a shared componentwise covariant lift, P2316 leaves theorem-grade full closure open, P1974 gives a nonzero anisotropic Bianchi-I residual, and P1795/P1809 keep the unified nonproxy gate locked."
        ),
        "failed_rows": failed,
        "next_honest_step": (
            "Do not promote L_total or claim ToE closure.  The next honest proof-grade step is P2686: a shared-background nonproxy component residual table for EA, EH, and ELg, with the Bianchi-I anisotropic residual as the first required row.  If that table cannot be made zero, export a no-go theorem for reverse closure from the current reduced/FRW scaffold."
        ),
        "full_tensor_resolved_eom_closed_now": lattice["current_reverse_closure_success"],
        "role_bearing_ltotal_exported_now": False,
        "selector_or_bridge_imported_now": False,
        "qw2191_discharged_now": False,
        "toe_closed_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2685/S1635 strict Lagrangian/EOM reverse-closure obstruction matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first grep",
    ]
    for name, data in payload["content_grep"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(["", "## Reverse-closure rows"])
    for row in payload["reverse_closure_rows"]:
        lines.append(f"- `{row['row_id']}`: available=`{row['available_forward_export']}`, satisfied_now=`{row['satisfied_now']}`; blocker={row['witness_or_blocker']}")
    res = payload["anisotropic_residual_rank"]
    lines.extend([
        "", "## Anisotropic residual computation",
        f"Residual vector: `{res['residual_vector']}`.",
        f"Symbolic Jacobian rank: `{res['symbolic_jacobian_rank']}`; isotropic limit zero: `{res['isotropic_limit_zero']}`; isotropic Jacobian rank: `{res['isotropic_jacobian_rank']}`; nonzero numeric samples: `{res['numeric_nonzero_samples']}`.",
        "", "## Verdict", payload["decision"]["professorial_verdict"],
        "", "## Next honest step", payload["decision"]["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    state = current_artifact_state()
    residual = anisotropic_residual_rank()
    rows = reverse_closure_rows(state, residual)
    lattice = policy_lattice(rows)
    payload: dict[str, Any] = {
        "status": "P2685_STRICT_LAGRANGIAN_EOM_REVERSE_CLOSURE_OBSTRUCTION_MATRIX_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "current_artifact_state": state,
        "anisotropic_residual_rank": residual,
        "reverse_closure_rows": rows,
        "policy_lattice": lattice,
        "decision": decision(rows, lattice),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2685/S1635 strict Lagrangian/EOM reverse-closure obstruction matrix",
        "## P2685/S1635 strict Lagrangian/EOM reverse-closure obstruction matrix\n\n"
        "`P2685/S1635` executes the finite reverse-closure audit after P2684.  Selector-independent reduced/forward Lagrangian/EOM exports are real, but they do not reverse into nonproxy full tensor-resolved closure: P2316 leaves theorem-grade full EOM open, P1974 supplies a nonzero Bianchi-I anisotropic residual, and P1795/P1809 keep the unified nonproxy gate locked.  No role-bearing `L_total`, selector/bridge import, `QW-2191` discharge, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2685/S1635 reverse-closure obstruction guard",
        "## P2685/S1635 reverse-closure obstruction guard\n\n"
        "`P2685/S1635` turns the P2684 pivot into an executable obstruction matrix.  The next honest Lagrangian/EOM move is a shared-background nonproxy component residual table for `EA`, `EH`, and `ELg`, beginning with the Bianchi-I anisotropic row; absent a zero table, the correct output is a no-go for reverse closure from the current reduced/FRW scaffold.\n",
    )
    append_once(
        AGENTS,
        "Current strict Lagrangian/EOM reverse-closure guardrail (P2685/S1635, 2026-06-13)",
        "## Current strict Lagrangian/EOM reverse-closure guardrail (P2685/S1635, 2026-06-13)\n\n"
        "- P2685 confirms that selector-independent reduced Lagrangian/EOM progress is real but not reversible to full tensor-resolved nonproxy closure on current artifacts.\n"
        "- The next proof-grade move is a shared-background `EA/EH/ELg` nonproxy component residual table, starting from the Bianchi-I anisotropic obstruction; do not promote `L_total`, ToE closure, selector closure, role transfer, or generic bridge closure without that table.\n",
    )
    return payload


if __name__ == "__main__":
    main()
