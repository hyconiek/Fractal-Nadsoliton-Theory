#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from itertools import product
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2606_s1556_strict_side_nonlinear_compression_residual_addition.json"
MD = GEN / "p2606_s1556_strict_side_nonlinear_compression_residual_addition.md"

SOURCE_FILES = {
    "P2603_FRACTAL_CODIMENSION_SLOPE_SOURCE": GEN / "p2603_s1553_nadsoliton_fractal_codimension_slope_source_theorem.json",
    "P2605_LINEAR_SLICE_COMPLETION_MAP": GEN / "p2605_s1555_legacy_strict_linear_slice_completion_map_evidence.json",
    "K1_KERNEL_SPLIT_NOTE": ROOT / "K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md",
}
BETA_TORS = 0.01
OMEGA = math.pi / 4.0
PHI = math.pi / 6.0
ETA_LINEAR = 1.0
ETA_NONLINEAR = 4.0 / 5.0
AUDIT_D_VALUES = list(range(1, 65))
REMAINING_AFTER_NONLINEAR = [
    "phase_topological_selector_data_certified",
    "strict_damping_role_transfer_theorem",
]
NEGATIVE_EXPORT_FLAGS = [
    "strict_side_residual_additions_exported",
    "phase_topological_selector_data_exported",
    "full_legacy_to_strict_completion_bridge_exported",
    "strict_damping_role_transfer_theorem_exported",
    "legacy_physical_role_transfer_exported",
    "role_bearing_ltotal_exported",
    "qw2191_discharged_by_this_addition",
    "toe_closure_claimed",
]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2606|S1556|nonlinear compression residual addition|strict-side nonlinear residual|eta four fifths bridge residual",
        "intended_research_nonduplication": "linear slice nonlinear residual|phase topological residual gate|post linear slice residual split|nonlinear strict-side addition",
        "precursor_chain": "P2603|S1553|P2605|S1555|D_f=9/5|eta=1|eta=4/5|K_strict_gate",
        "guardrails": "role-transfer theorem|role-bearing L_total|QW-2191|ToE closure|kernel-split guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def strict_kernel(d: int, eta: float) -> float:
    return math.cos(OMEGA * d + PHI) / (1.0 + BETA_TORS * (d ** eta))


def nonlinear_residual_audit() -> dict[str, Any]:
    rows = []
    linear_values = []
    nonlinear_values = []
    residual_values = []
    ratio_values = []
    for d in AUDIT_D_VALUES:
        linear = strict_kernel(d, ETA_LINEAR)
        nonlinear = strict_kernel(d, ETA_NONLINEAR)
        residual = nonlinear - linear
        linear_values.append(linear)
        nonlinear_values.append(nonlinear)
        residual_values.append(residual)
        if abs(linear) > 1e-12:
            ratio_values.append(nonlinear / linear)
        rows.append({
            "d": d,
            "strict_linear_slice_eta_1": linear,
            "strict_nonlinear_eta_4_over_5": nonlinear,
            "nonlinear_residual_eta_4_over_5_minus_eta_1": residual,
            "abs_residual": abs(residual),
        })
    dot_ln = sum(l * n for l, n in zip(linear_values, nonlinear_values))
    dot_ll = sum(l * l for l in linear_values)
    best_scalar = dot_ln / dot_ll
    scalar_residuals = [n - best_scalar * l for l, n in zip(linear_values, nonlinear_values)]
    return {
        "residual_formula": "R_nonlinear(d)=cos(omega*d+phi)/(1+beta_tors*d^(4/5))-cos(omega*d+phi)/(1+beta_tors*d)",
        "eta_linear_slice": ETA_LINEAR,
        "eta_nonlinear_codimension": ETA_NONLINEAR,
        "beta_value": BETA_TORS,
        "audit_d_values": AUDIT_D_VALUES,
        "pointwise_rows": rows,
        "max_abs_nonlinear_residual": max(abs(x) for x in residual_values),
        "nonzero_residual_count": sum(1 for x in residual_values if abs(x) > 1e-12),
        "best_scalar_fit_nonlinear_as_scaled_linear": best_scalar,
        "max_abs_residual_after_best_scalar_fit": max(abs(x) for x in scalar_residuals),
        "ratio_spread_nonlinear_over_linear": max(ratio_values) - min(ratio_values),
        "nonlinear_residual_not_absorbed_by_amplitude_rescaling": max(abs(x) for x in scalar_residuals) > 1e-4,
        "component_interpretation": "The eta=4/5 codimension law supplies a genuine nonlinear damping/compression denominator residual beyond the eta=1 legacy linear slice; it is not a full bridge because phase/topological selector data remain separate.",
    }


def remaining_rows_after_nonlinear() -> list[dict[str, Any]]:
    rows = []
    for values in product([False, True], repeat=len(REMAINING_AFTER_NONLINEAR)):
        assignment = dict(zip(REMAINING_AFTER_NONLINEAR, values))
        assignment["linear_slice_completion_map_evidence"] = True
        assignment["nonlinear_compression_residual_addition"] = True
        role_bearing = all(assignment.values())
        missing = [key for key, value in assignment.items() if not value]
        rows.append({
            "assignment": assignment,
            "role_bearing_ltotal_accepts_after_nonlinear_component": role_bearing,
            "missing_remaining_gates": missing,
            "missing_remaining_gate_count": len(missing),
        })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2603_payload = load_json(SOURCE_FILES["P2603_FRACTAL_CODIMENSION_SLOPE_SOURCE"])
    p2605_payload = load_json(SOURCE_FILES["P2605_LINEAR_SLICE_COMPLETION_MAP"])
    p2603 = theorem(p2603_payload, "nadsoliton_fractal_codimension_slope_source_theorem")
    p2605 = theorem(p2605_payload, "legacy_strict_linear_slice_completion_map_evidence")
    audit = nonlinear_residual_audit()
    rows = remaining_rows_after_nonlinear()
    accepting = [row for row in rows if row["role_bearing_ltotal_accepts_after_nonlinear_component"]]
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2606_T1_strict_side_nonlinear_compression_residual_addition",
        "audited_chain": ["P2603/S1553", "P2605/S1555", "K1"],
        "nonlinear_residual_statement": (
            "Using the codimension exponent eta=4/5 sourced in P2603, the strict gate denominator differs from the eta=1 legacy linear slice by a nonzero nonlinear compression residual R_nonlinear(d). This supplies one strict-side residual component but not the phase/topological selector component or role transfer."
        ),
        "p2603_codimension_slope_inherited": p2603.get("slope_value_or_prime_anchor_source_exported") is True,
        "p2605_linear_slice_evidence_inherited": p2605.get("linear_slice_completion_map_evidence_exported") is True,
        "nonlinear_compression_residual_addition_exported": True,
        "strict_side_nonlinear_residual_audit": audit,
        "post_nonlinear_remaining_bridge_matrix": {
            "remaining_gates_after_nonlinear_component": REMAINING_AFTER_NONLINEAR,
            "truth_table_rows": rows,
            "truth_table_row_count": len(rows),
            "accepting_row_count": len(accepting),
            "current_assignment_after_p2606": {
                "linear_slice_completion_map_evidence": True,
                "nonlinear_compression_residual_addition": True,
                "phase_topological_selector_data_certified": False,
                "strict_damping_role_transfer_theorem": False,
            },
            "current_role_bearing_ltotal_accepts": False,
            "remaining_missing_gate_count_after_p2606": 2,
        },
        "recommended_next_honest_step": (
            "Next certify the remaining strict-side phase/topological selector data, not APD/moment/Sturm and not legacy physical-role transfer. Role transfer remains inadmissible until that selector data is supplied."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2603_codimension_inherited": theorem_export["p2603_codimension_slope_inherited"],
        "p2605_linear_slice_inherited": theorem_export["p2605_linear_slice_evidence_inherited"],
        "nonlinear_component_exported": theorem_export["nonlinear_compression_residual_addition_exported"],
        "nonlinear_residual_nonzero": audit["nonzero_residual_count"] > 0 and audit["max_abs_nonlinear_residual"] > 1e-4,
        "nonlinear_residual_not_scalar_rescaling": audit["nonlinear_residual_not_absorbed_by_amplitude_rescaling"],
        "remaining_truth_table_has_four_rows": theorem_export["post_nonlinear_remaining_bridge_matrix"]["truth_table_row_count"] == 4,
        "remaining_accepting_row_count_one": theorem_export["post_nonlinear_remaining_bridge_matrix"]["accepting_row_count"] == 1,
        "current_role_bearing_ltotal_blocked": theorem_export["post_nonlinear_remaining_bridge_matrix"]["current_role_bearing_ltotal_accepts"] is False,
        "two_remaining_gates_after_p2606": theorem_export["post_nonlinear_remaining_bridge_matrix"]["remaining_missing_gate_count_after_p2606"] == 2,
        "no_full_strict_side_residuals_exported": theorem_export["strict_side_residual_additions_exported"] is False,
        "no_phase_topological_exported": theorem_export["phase_topological_selector_data_exported"] is False,
        "no_role_transfer_exported": theorem_export["strict_damping_role_transfer_theorem_exported"] is False,
        "no_role_bearing_ltotal_exported": theorem_export["role_bearing_ltotal_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_addition"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2606",
        "stage_id": "S1556",
        "status": "P2606_NONLINEAR_COMPRESSION_RESIDUAL_COMPONENT_EXPORTED_PHASE_TOPOLOGICAL_AND_ROLE_TRANSFER_STILL_BLOCKED_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_side_nonlinear_compression_residual_addition": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2603_FRACTAL_CODIMENSION_SLOPE_SOURCE": sha256_json(p2603_payload),
                "P2605_LINEAR_SLICE_COMPLETION_MAP": sha256_json(p2605_payload),
                "K1_KERNEL_SPLIT_NOTE": sha256_bytes(Path(SOURCE_FILES["K1_KERNEL_SPLIT_NOTE"]).read_bytes()),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_side_nonlinear_compression_residual_addition"]["theorem_export"]
    audit = t["strict_side_nonlinear_residual_audit"]
    remaining = t["post_nonlinear_remaining_bridge_matrix"]
    lines = [
        "# P2606/S1556 strict-side nonlinear compression residual addition", "",
        f"Status: `{payload['status']}`", "", "## Residual statement", "",
        t["nonlinear_residual_statement"], "", "## Computed consequences", "",
        f"- Nonlinear residual component exported: `{t['nonlinear_compression_residual_addition_exported']}`.",
        f"- Max nonlinear residual: `{audit['max_abs_nonlinear_residual']}`.",
        f"- Residual count above tolerance: `{audit['nonzero_residual_count']}`.",
        f"- Max residual after best scalar fit: `{audit['max_abs_residual_after_best_scalar_fit']}`.",
        f"- Remaining gates: `{remaining['remaining_gates_after_nonlinear_component']}`.",
        f"- Current role-bearing L_total accepts: `{remaining['current_role_bearing_ltotal_accepts']}`.", "",
        "## Scope guards", "",
        "This exports only the nonlinear compression residual component. It does not export the full strict-side residual additions package, phase/topological selector data, role-transfer theorem, role-bearing `L_total`, QW-2191 discharge, or ToE closure.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Fingerprint", "",
        f"`{payload['strict_side_nonlinear_compression_residual_addition']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2606/S1556 strict-side nonlinear compression residual addition

`P2606/S1556` supplies one strict-side residual component beyond the P2605 linear slice: with the P2603 codimension exponent `eta=4/5`, the strict denominator differs nontrivially from the legacy `eta=1` denominator and the residual is not absorbed by a scalar amplitude rescaling.  This still does not certify the full strict-side residual additions package, because phase/topological selector data and role-transfer remain separate gates.
""".strip()
    lag_section = """
## P2606/S1556 nonlinear residual Ltotal guard

`P2606/S1556` keeps the damping/compression term non-role-bearing in `L_total`: the nonlinear compression residual component is now isolated, but phase/topological selector data and a strict damping role-transfer theorem remain required before any role-bearing promotion.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2606/S1556 strict-side nonlinear compression residual addition", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2606/S1556 nonlinear residual Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
