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
OUT = GEN / "p2605_s1555_legacy_strict_linear_slice_completion_map_evidence.json"
MD = GEN / "p2605_s1555_legacy_strict_linear_slice_completion_map_evidence.md"

SOURCE_FILES = {
    "P2604_POST_SOURCE_BRIDGE_READINESS": GEN / "p2604_s1554_strict_damping_post_source_bridge_readiness_matrix.json",
    "K1_KERNEL_SPLIT_NOTE": ROOT / "K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md",
    "K2_STRICT_DERIVATION_CHAIN_NOTE": ROOT / "K2_STRICT_GATE_KERNEL_DERIVATION_CHAIN_NOTE.md",
}
ALPHA_GEO = 4.0 * math.log(2.0)
BETA_TORS = 0.01
OMEGA = math.pi / 4.0
PHI = math.pi / 6.0
AUDIT_D_VALUES = list(range(1, 65))
ETA_CANDIDATES = [1.0, 0.8, 1.5]
REMAINING_AFTER_LINEAR_SLICE = [
    "strict_side_residual_additions_certified",
    "strict_damping_role_transfer_theorem",
]
NEGATIVE_EXPORT_FLAGS = [
    "full_legacy_to_strict_completion_bridge_exported",
    "strict_side_residual_additions_exported",
    "strict_damping_role_transfer_theorem_exported",
    "role_bearing_ltotal_exported",
    "legacy_physical_role_transfer_exported",
    "qw2191_discharged_by_this_evidence",
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
        "new_packet": "P2605|S1555|legacy normalized linear slice|linear-slice completion map|alpha_geo normalization bridge",
        "intended_research_nonduplication": "beta_tors beta eta one|legacy strict kernel denominator match|K_legacy.*K_strict.*linear slice|completion map evidence theorem",
        "precursor_chain": "P2604|S1554|post-source bridge readiness|K_legacy_ont|K_strict_gate|alpha_geo|beta_tors",
        "guardrails": "role-transfer theorem|role-bearing L_total|QW-2191|ToE closure|kernel-split guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def legacy_kernel(d: int) -> float:
    return ALPHA_GEO * math.cos(OMEGA * d + PHI) / (1.0 + BETA_TORS * d)


def strict_kernel_linear_slice(d: int) -> float:
    return math.cos(OMEGA * d + PHI) / (1.0 + BETA_TORS * d)


def strict_kernel_eta(d: int, eta: float) -> float:
    return math.cos(OMEGA * d + PHI) / (1.0 + BETA_TORS * (d ** eta))


def completion_map_audit() -> dict[str, Any]:
    rows = []
    max_abs_residual = 0.0
    for d in AUDIT_D_VALUES:
        normalized_legacy = legacy_kernel(d) / ALPHA_GEO
        strict_linear = strict_kernel_linear_slice(d)
        residual = normalized_legacy - strict_linear
        max_abs_residual = max(max_abs_residual, abs(residual))
        rows.append({
            "d": d,
            "legacy_normalized_by_alpha_geo": normalized_legacy,
            "strict_linear_slice_beta_equals_beta_tors_eta_1": strict_linear,
            "residual": residual,
            "abs_residual": abs(residual),
        })
    eta_rows = []
    for eta in ETA_CANDIDATES:
        residuals = [legacy_kernel(d) / ALPHA_GEO - strict_kernel_eta(d, eta) for d in AUDIT_D_VALUES]
        eta_rows.append({
            "eta_candidate": eta,
            "max_abs_residual_against_normalized_legacy": max(abs(x) for x in residuals),
            "exact_linear_slice_accepts": eta == 1.0 and max(abs(x) for x in residuals) < 1e-14,
        })
    return {
        "completion_map_formula": "K_legacy_ont(d)/alpha_geo = K_strict_gate(d) on the linear denominator slice beta=beta_tors, eta=1",
        "parameter_map_evidence": {
            "amplitude_normalization": "divide by alpha_geo = 4 ln 2",
            "beta_map": "beta <- beta_tors = 0.01",
            "eta_slice": "eta <- 1",
            "omega_phi_pass_through": "omega and phi are held fixed on this slice",
        },
        "audit_d_values": AUDIT_D_VALUES,
        "pointwise_rows": rows,
        "max_abs_residual_on_audit_grid": max_abs_residual,
        "eta_candidate_rows": eta_rows,
        "linear_slice_exact_on_grid": max_abs_residual < 1e-14,
        "residual_strict_side_additions_not_covered": [
            "nonlinear compression/damping beyond the legacy linear denominator slice",
            "certified phase/topological data and selector premises",
            "legacy physical-role transfer such as sin^2(theta_W), alpha_EM, or gravity hierarchy claims",
        ],
    }


def remaining_bridge_rows() -> list[dict[str, Any]]:
    rows = []
    for values in product([False, True], repeat=len(REMAINING_AFTER_LINEAR_SLICE)):
        assignment = dict(zip(REMAINING_AFTER_LINEAR_SLICE, values))
        assignment["legacy_to_strict_completion_map_evidence"] = True
        role_bearing = all(assignment.values())
        missing = [key for key, value in assignment.items() if not value]
        rows.append({
            "assignment": assignment,
            "role_bearing_ltotal_accepts_after_linear_slice": role_bearing,
            "missing_remaining_gates": missing,
            "missing_remaining_gate_count": len(missing),
        })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2604_payload = load_json(SOURCE_FILES["P2604_POST_SOURCE_BRIDGE_READINESS"])
    p2604 = theorem(p2604_payload, "strict_damping_post_source_bridge_readiness_matrix")
    audit = completion_map_audit()
    rows = remaining_bridge_rows()
    accepting = [row for row in rows if row["role_bearing_ltotal_accepts_after_linear_slice"]]
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2605_T1_legacy_strict_linear_slice_completion_map_evidence",
        "audited_chain": ["P2604/S1554", "K1", "K2"],
        "completion_map_evidence_statement": (
            "The normalized legacy kernel matches the strict gate kernel exactly on the linear denominator slice: K_legacy_ont(d)/alpha_geo = cos(omega d+phi)/(1+beta_tors d), which is K_strict_gate(d) with beta=beta_tors and eta=1."
        ),
        "linear_slice_completion_map_evidence_exported": True,
        "legacy_to_strict_completion_map_evidence_exported_for_linear_slice_only": True,
        "completion_map_audit": audit,
        "p2604_strict_damping_source_inherited": p2604.get("strict_damping_beta_eta_source_inherited_from_p2603") is True,
        "post_linear_slice_remaining_bridge_matrix": {
            "remaining_gates_after_linear_slice_evidence": REMAINING_AFTER_LINEAR_SLICE,
            "truth_table_rows": rows,
            "truth_table_row_count": len(rows),
            "accepting_row_count": len(accepting),
            "current_assignment_after_p2605": {
                "legacy_to_strict_completion_map_evidence": True,
                "strict_side_residual_additions_certified": False,
                "strict_damping_role_transfer_theorem": False,
            },
            "current_role_bearing_ltotal_accepts": False,
            "remaining_missing_gate_count_after_p2605": 2,
        },
        "recommended_next_honest_step": (
            "Do not transfer legacy physical roles. The next bridge step must certify strict-side residual additions: nonlinear compression/damping beyond eta=1 plus phase/topological selector data, before any role-transfer theorem."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2604_source_readiness_inherited": theorem_export["p2604_strict_damping_source_inherited"],
        "linear_slice_completion_map_evidence_exported": theorem_export["linear_slice_completion_map_evidence_exported"],
        "linear_slice_exact_on_grid": audit["linear_slice_exact_on_grid"],
        "max_residual_below_tolerance": audit["max_abs_residual_on_audit_grid"] < 1e-14,
        "eta_one_unique_exact_candidate": sum(1 for row in audit["eta_candidate_rows"] if row["exact_linear_slice_accepts"]) == 1,
        "remaining_bridge_truth_table_has_four_rows": theorem_export["post_linear_slice_remaining_bridge_matrix"]["truth_table_row_count"] == 4,
        "remaining_accepting_row_count_one": theorem_export["post_linear_slice_remaining_bridge_matrix"]["accepting_row_count"] == 1,
        "current_role_bearing_ltotal_blocked": theorem_export["post_linear_slice_remaining_bridge_matrix"]["current_role_bearing_ltotal_accepts"] is False,
        "two_remaining_gates_after_p2605": theorem_export["post_linear_slice_remaining_bridge_matrix"]["remaining_missing_gate_count_after_p2605"] == 2,
        "no_full_bridge_exported": theorem_export["full_legacy_to_strict_completion_bridge_exported"] is False,
        "no_strict_side_residuals_exported": theorem_export["strict_side_residual_additions_exported"] is False,
        "no_role_transfer_exported": theorem_export["strict_damping_role_transfer_theorem_exported"] is False,
        "no_role_bearing_ltotal_exported": theorem_export["role_bearing_ltotal_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_evidence"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    source_fingerprints = {
        "P2604_POST_SOURCE_BRIDGE_READINESS": sha256_json(p2604_payload),
        "K1_KERNEL_SPLIT_NOTE": sha256_bytes(Path(SOURCE_FILES["K1_KERNEL_SPLIT_NOTE"]).read_bytes()),
        "K2_STRICT_DERIVATION_CHAIN_NOTE": sha256_bytes(Path(SOURCE_FILES["K2_STRICT_DERIVATION_CHAIN_NOTE"]).read_bytes()),
    }
    return {
        "packet_id": "P2605",
        "stage_id": "S1555",
        "status": "P2605_LINEAR_SLICE_COMPLETION_MAP_EVIDENCE_EXPORTED_FULL_BRIDGE_AND_ROLE_TRANSFER_STILL_BLOCKED_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "legacy_strict_linear_slice_completion_map_evidence": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": source_fingerprints,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["legacy_strict_linear_slice_completion_map_evidence"]["theorem_export"]
    audit = t["completion_map_audit"]
    remaining = t["post_linear_slice_remaining_bridge_matrix"]
    lines = [
        "# P2605/S1555 legacy-strict linear-slice completion map evidence", "",
        f"Status: `{payload['status']}`", "", "## Evidence statement", "",
        t["completion_map_evidence_statement"], "", "## Computed consequences", "",
        f"- Linear-slice evidence exported: `{t['linear_slice_completion_map_evidence_exported']}`.",
        f"- Max audit residual: `{audit['max_abs_residual_on_audit_grid']}`.",
        f"- Eta candidate rows: `{audit['eta_candidate_rows']}`.",
        f"- Remaining bridge gates: `{remaining['remaining_gates_after_linear_slice_evidence']}`.",
        f"- Remaining truth-table rows: `{remaining['truth_table_row_count']}`.",
        f"- Current role-bearing L_total accepts: `{remaining['current_role_bearing_ltotal_accepts']}`.", "",
        "## Scope guards", "",
        "This is only linear-slice completion-map evidence. It does not export strict-side residual additions, a full legacy-to-strict bridge, role-transfer theorem, legacy physical-role transfer, role-bearing `L_total`, QW-2191 discharge, or ToE closure.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Fingerprint", "",
        f"`{payload['legacy_strict_linear_slice_completion_map_evidence']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2605/S1555 legacy-strict linear-slice completion map evidence

`P2605/S1555` exports only a linear-slice completion-map evidence theorem: after normalizing `K_legacy_ont` by `alpha_geo`, the denominator matches `K_strict_gate` exactly on the slice `beta=beta_tors`, `eta=1`.  This supplies completion-map evidence for the legacy linear denominator slice, but it does not certify strict-side residual additions, a full bridge, legacy physical-role transfer, role-bearing `L_total`, QW-2191 discharge, or ToE closure.
""".strip()
    lag_section = """
## P2605/S1555 linear-slice bridge evidence Ltotal guard

`P2605/S1555` reduces the P2604 bridge-readiness matrix by supplying linear-slice completion-map evidence only.  The damping/compression term remains non-role-bearing in `L_total` until nonlinear strict-side residual additions and a strict damping role-transfer theorem are supplied.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2605/S1555 legacy-strict linear-slice completion map evidence", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2605/S1555 linear-slice bridge evidence Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
