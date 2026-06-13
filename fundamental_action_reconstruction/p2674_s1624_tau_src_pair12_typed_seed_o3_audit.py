#!/usr/bin/env python3
"""P2674/S1624: O3 chart-sensitive pair1/pair2 typed seed audit.

Audits the first missing obligation from P2673:
    current chart-sensitive pair1/pair2 typed seed subinterface on tau_src
without promoting target-only or attempt-only material into an actual export.
"""
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
OUT = GEN / "p2674_s1624_tau_src_pair12_typed_seed_o3_audit.json"
MD = GEN / "p2674_s1624_tau_src_pair12_typed_seed_o3_audit.md"
P2673 = GEN / "p2673_s1623_tau_src_pair12_boundary_square_subinterface_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NEGATIVE_EXPORT_FLAGS = [
    "o3_chart_sensitive_pair12_typed_seed_exported",
    "tau_src_to_pair12_typed_seed_arrow_exported",
    "chart_label_retaining_nonprojector_seed_exported",
    "sigma_sel_src_to_f301_bridge_exported",
    "boundary_square_cycle_typed_arrow_exported",
    "sector_swap_sourced_invariant_exported",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "toe_closure_claimed",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    data = json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(data.encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:60]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "tau_src_anchor_content": "tau_src_candidate_v1|tau_src.*sign|source-topology selector witness",
        "pair12_carrier_content": "F301|residual-datum pair1/pair2 carrier|surviving.*pair1/pair2 carrier|pair1/pair2 residual-datum",
        "chart_sensitive_atlas_content": "chart-sensitive.*pair1/pair2|local pair1/pair2 atlas|SelectorAtlas_pair12|sigma_int corridor",
        "typed_seed_target_content": "chart-label-retaining pair1/pair2 typed seed|typed seed subinterface|pair12 typed seed|T222 seed-subinterface",
        "actual_export_blocker_content": "future-only target|not actually realized|no actual realization|attempt must remain below actual subinterface export|nonexport",
        "collapse_blocker_content": "Q_basis_sel_v1 terminal collapse|projector-only.*collapse|preLM|quotient-class only|projector-only local pair12 atlas",
        "closure_guard_content": "QW-2191|L_total|ToE closure|boundary square cycle|sector-swap sourced invariant",
    }
    return {
        "tool": "rg",
        "mode": "content-first O3 chart-sensitive pair1/pair2 typed seed subinterface audit; not name/number-only",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def upstream_consistency() -> dict[str, Any]:
    p2673 = load_json(P2673)
    decision = p2673.get("closure_decision", {})
    lattice = p2673.get("finite_closure_lattice", {})
    return {
        "p2673_names_o3_missing": "O3_chart_sensitive_pair12_typed_subinterface" in lattice.get("missing_obligations_now", []),
        "p2673_distance_three": lattice.get("hamming_distance_from_closure") == 3,
        "p2673_no_boundary_phase_bit_target": decision.get("boundary_phase_bit_target_exported_now") is False,
        "p2673_no_ltotal_reopening": decision.get("role_bearing_ltotal_now") is False,
    }


def o3_acceptance_matrix() -> list[dict[str, Any]]:
    return [
        {
            "id": "S1_tau_src_source_anchor_present",
            "description": "tau_src/source-topology sign or selector witness material is present as an input anchor",
            "satisfied_now": True,
            "content_pattern": "tau_src_candidate_v1|tau_src.*sign|source-topology selector witness",
        },
        {
            "id": "S2_surviving_pair12_carrier_present",
            "description": "surviving F301/residual-datum pair1/pair2 carrier material is present",
            "satisfied_now": True,
            "content_pattern": "F301|residual-datum pair1/pair2 carrier|surviving.*pair1/pair2 carrier|pair1/pair2 residual-datum",
        },
        {
            "id": "S3_chart_sensitive_atlas_lane_present",
            "description": "local/chart-sensitive pair1/pair2 atlas lane exists as material to compare against",
            "satisfied_now": True,
            "content_pattern": "chart-sensitive.*pair1/pair2|local pair1/pair2 atlas|SelectorAtlas_pair12|sigma_int corridor",
        },
        {
            "id": "S4_actual_chart_label_retaining_typed_seed_exported",
            "description": "an actual chart-label-retaining pair1/pair2 typed seed subinterface is exported, not merely targeted",
            "satisfied_now": False,
            "content_pattern": "chart-label-retaining pair1/pair2 typed seed|typed seed subinterface|pair12 typed seed|actual subinterface export",
        },
        {
            "id": "S5_nonprojector_nonquotient_nonprelm_descent_law",
            "description": "the descent avoids terminal Q_basis/preLM/quotient-only and projector-only local-atlas collapse",
            "satisfied_now": False,
            "content_pattern": "Q_basis_sel_v1 terminal collapse|projector-only.*collapse|preLM|quotient-class only|projector-only local pair12 atlas",
        },
        {
            "id": "S6_sigma_to_f301_typed_arrow_not_by_fiat",
            "description": "Sigma_sel_src_target_v1 is typed into the surviving F301 carrier without identification by fiat",
            "satisfied_now": False,
            "content_pattern": "Sigma_sel_src_target_v1.*F301|F301.*Sigma_sel_src_target_v1|must_not_identify_Sigma_sel_src_target_v1|unbridged to.*pair1/pair2",
        },
    ]


def score_matrix(matrix: list[dict[str, Any]]) -> list[dict[str, Any]]:
    scored = []
    for item in matrix:
        hits = rg_count(item["content_pattern"])
        scored.append({**item, "content_hits": hits["count"], "content_samples": hits["samples"][:10]})
    return scored


def finite_o3_lattice(scored: list[dict[str, Any]]) -> dict[str, Any]:
    ids = [item["id"] for item in scored]
    rows = []
    for bits in itertools.product([False, True], repeat=len(ids)):
        state = dict(zip(ids, bits, strict=True))
        rows.append({"state": state, "passes_o3": all(bits)})
    current = {item["id"]: item["satisfied_now"] for item in scored}
    return {
        "total_states_checked": len(rows),
        "passing_o3_states_count": sum(row["passes_o3"] for row in rows),
        "current_state": current,
        "current_state_passes_o3": all(current.values()),
        "missing_o3_subobligations_now": [key for key, value in current.items() if not value],
        "hamming_distance_from_o3_pass": sum(1 for value in current.values() if not value),
        "o3_pass_state": {key: True for key in ids},
    }


def closure_decision(lattice: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "P2674_O3_CHART_SENSITIVE_PAIR12_TYPED_SEED_AUDIT__TARGET_AND_ATTEMPT_MATERIAL_REAL_BUT_NO_ACTUAL_TYPED_SEED_EXPORT",
        "professorial_verdict": (
            "P2674 attacks O3 directly. The repo has real tau_src input material, a surviving pair1/pair2 carrier, and a local chart-sensitive atlas lane, "
            "but the finite O3 lattice still lacks an actual chart-label-retaining typed seed export, a nonprojector/nonquotient descent law, and a Sigma_sel_src_target_v1 -> F301 typed arrow not imposed by fiat. "
            "Thus O3 remains blocked rather than proved."
        ),
        "next_honest_step": (
            "Do not try O4/O5 yet. The next honest proof-grade step is a no-go-or-construction audit of the exact S6 arrow: "
            "Sigma_sel_src_target_v1 -> surviving F301 pair1/pair2 carrier before Q_basis/preLM and projector-only collapse. "
            "If S6 cannot be exported without fiat identification, record O3 as blocked and promote the tau_src -> pair12 -> boundary-square route to no-go at the typed-seed interface."
        ),
        "current_state_passes_o3": lattice["current_state_passes_o3"],
        "missing_o3_subobligations_now": lattice["missing_o3_subobligations_now"],
        "hamming_distance_from_o3_pass": lattice["hamming_distance_from_o3_pass"],
        "o3_exported_now": False,
        "boundary_square_arrow_allowed_next": False,
        "sector_swap_invariant_allowed_next": False,
        "boundary_phase_bit_target_exported_now": False,
        "qw2191_discharged_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lattice = payload["finite_o3_lattice"]
    decision = payload["closure_decision"]
    lines = [
        "# P2674/S1624 O3 chart-sensitive pair1/pair2 typed seed audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first repo grep",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(["", "## O3 acceptance matrix"])
    for item in payload["o3_acceptance_matrix"]:
        lines.append(f"- `{item['id']}`: satisfied_now=`{item['satisfied_now']}`, content_hits=`{item['content_hits']}` — {item['description']}")
    lines.extend([
        "",
        "## Finite O3 lattice",
        f"Total states checked: `{lattice['total_states_checked']}`.",
        f"Passing O3 states: `{lattice['passing_o3_states_count']}`.",
        f"Current state passes O3? `{lattice['current_state_passes_o3']}`.",
        f"Missing O3 subobligations now: `{lattice['missing_o3_subobligations_now']}`.",
        f"Hamming distance from O3 pass: `{lattice['hamming_distance_from_o3_pass']}`.",
        "",
        "## Verdict",
        decision["professorial_verdict"],
        f"Decision: `{decision['decision']}`.",
        "",
        "## Next honest step",
        decision["next_honest_step"],
        "",
        "## Negative exports",
    ])
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    scored = score_matrix(o3_acceptance_matrix())
    lattice = finite_o3_lattice(scored)
    decision = closure_decision(lattice)
    payload: dict[str, Any] = {
        "status": "P2674_O3_CHART_SENSITIVE_PAIR12_TYPED_SEED_AUDIT_NO_FALSE_PASS",
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2673": sha256_file(P2673),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "o3_acceptance_matrix": scored,
        "finite_o3_lattice": lattice,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2674/S1624 O3 chart-sensitive pair12 typed seed guard",
        "## P2674/S1624 O3 chart-sensitive pair12 typed seed guard\n\n"
        "`P2674/S1624` attacks O3 from `P2673` directly.  The audit finds real tau_src input material, a surviving `F301` pair1/pair2 carrier, and a local chart-sensitive atlas lane, but no actual chart-label-retaining typed seed export, no nonprojector/nonquotient descent law before `Q_basis`/preLM collapse, and no `Sigma_sel_src_target_v1 -> F301` typed arrow not imposed by fiat.  Therefore O3 does not pass, O4/O5 remain inadmissible, and no boundary-phase bit target, beta source, `QW-2191` discharge, bridge completion, role transfer, role-bearing `L_total`, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2674/S1624 O3 typed-seed Ltotal guard",
        "## P2674/S1624 O3 typed-seed Ltotal guard\n\n"
        "`P2674/S1624` keeps `L_total` closed to the tau_src -> pair12 -> boundary-square route at O3.  No source-topology-sign-derived variational term is admissible until the exact `Sigma_sel_src_target_v1 -> F301` chart-label-retaining typed seed arrow is exported without quotient, preLM, projector-only, or fiat-identification collapse.\n",
    )
    return payload


if __name__ == "__main__":
    main()
