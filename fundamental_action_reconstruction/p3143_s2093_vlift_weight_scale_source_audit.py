"""P3143/S2093: V_lift weight/scale strict-source audit.

P3142 left one narrow honest continuation: try exactly one source audit for the
weights/scale of the axiom-branch selector functional.  This packet constructs a
finite candidate matrix from repo-backed scalar lanes and checks whether any
candidate can source mu, w_theta, w_s, or kappa without importing selector
axioms, local charts, unit conventions, or closed replay lanes.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3143_s2093_vlift_weight_scale_source_audit.json"
MD = GEN / "p3143_s2093_vlift_weight_scale_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P3142": GEN / "p3142_s2092_axiom_selector_field_variational_lift.json",
    "P3141": GEN / "p3141_s2091_axiom_selector_potential_downstream_readiness.json",
    "P2689": GEN / "p2689_s1639_entropy_reference_cell_bit_to_length_uv_unit_obligation_matrix.json",
    "P2691": GEN / "p2691_s1641_alpha_geo_role_safe_amplitude_source_audit.json",
    "P2692": GEN / "p2692_s1642_target_independent_positive_beta_zbeta_source_audit.json",
    "P3139": GEN / "p3139_s2089_dhl_lane_no_new_frontier_reconciliation.json",
}

TARGETS = ["mu", "w_theta", "w_s", "kappa"]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def candidate_sources() -> list[dict[str, Any]]:
    alpha_geo = 4.0 * math.log(2.0)
    return [
        {
            "id": "C1_entropy_bit_cell",
            "formula": "one-bit entropy reference cell -> positive scale",
            "value": 1.0,
            "positive_scalar": True,
            "repo_backing": "P2689",
            "blocked_by": ["intrinsic entropy reference missing", "selector-free bit target missing", "bit-to-length/action map missing", "scale-orbit quotient missing"],
            "imports_closed_lane": True,
            "unit_bearing": False,
            "global_Z12_compatible": True,
            "strict_source_law": False,
        },
        {
            "id": "C2_alpha_geo_shape_norm",
            "formula": "alpha_geo=4 ln 2 scalar-shape normalization",
            "value": alpha_geo,
            "positive_scalar": alpha_geo > 0,
            "repo_backing": "P2691",
            "blocked_by": ["role-safe amplitude absorption source missing", "physical-role safety theorem missing", "APD/Lagrangian dynamical source missing"],
            "imports_closed_lane": True,
            "unit_bearing": False,
            "global_Z12_compatible": True,
            "strict_source_law": False,
        },
        {
            "id": "C3_beta_z_beta_positive_orbit",
            "formula": "positive beta/Z_beta orbit representative",
            "value": 1.0,
            "positive_scalar": True,
            "repo_backing": "P2692",
            "blocked_by": ["target-independent positive beta/Z_beta theorem missing", "tail-ratio inversion is target-dependent"],
            "imports_closed_lane": True,
            "unit_bearing": False,
            "global_Z12_compatible": True,
            "strict_source_law": False,
        },
        {
            "id": "C4_vlift_hessian_self_normalization",
            "formula": "normalize weights from P3142 Hessian eigenvalue ratios",
            "value": 2.0,
            "positive_scalar": True,
            "repo_backing": "P3142",
            "blocked_by": ["circular: uses assumed V_lift weights to source V_lift weights", "local chart at r0 imported", "no independent unit measure"],
            "imports_closed_lane": False,
            "unit_bearing": False,
            "global_Z12_compatible": False,
            "strict_source_law": False,
        },
        {
            "id": "C5_dhl_receiver_obligation_gap",
            "formula": "D_HL receiver-to-source obstruction as negative source datum",
            "value": 0.0,
            "positive_scalar": False,
            "repo_backing": "P3139",
            "blocked_by": ["absolute support origin missing", "unpaired lambda missing", "import-free source law missing", "variational/unit coupling missing"],
            "imports_closed_lane": False,
            "unit_bearing": False,
            "global_Z12_compatible": False,
            "strict_source_law": False,
        },
    ]


def row_for(candidate: dict[str, Any], target: str) -> dict[str, Any]:
    # A candidate must satisfy all obligations to source one V_lift coefficient.
    obligations = {
        "positive_value": bool(candidate["positive_scalar"]),
        "strict_source_law": bool(candidate["strict_source_law"]),
        "unit_bearing": bool(candidate["unit_bearing"]),
        "global_Z12_compatible": bool(candidate["global_Z12_compatible"]),
        "noncircular_for_target": not (candidate["id"] == "C4_vlift_hessian_self_normalization" and target in TARGETS),
        "no_closed_lane_replay": not bool(candidate["imports_closed_lane"]),
    }
    accepted = all(obligations.values())
    return {
        "candidate_id": candidate["id"],
        "target": target,
        "repo_backing": candidate["repo_backing"],
        "value": candidate["value"],
        "obligations": obligations,
        "accepted_as_strict_source": accepted,
        "blocking_reasons": [] if accepted else candidate["blocked_by"],
    }


def build_payload() -> dict[str, Any]:
    inputs = {name: load_json(path) for name, path in INPUTS.items()}
    candidates = candidate_sources()
    rows = [row_for(c, target) for c in candidates for target in TARGETS]
    accepted = [r for r in rows if r["accepted_as_strict_source"]]
    source_defect_table = [
        {"obligation": "positive scalar value", "rows_passing": sum(r["obligations"]["positive_value"] for r in rows), "rows_total": len(rows)},
        {"obligation": "strict source law", "rows_passing": sum(r["obligations"]["strict_source_law"] for r in rows), "rows_total": len(rows)},
        {"obligation": "unit-bearing normalization", "rows_passing": sum(r["obligations"]["unit_bearing"] for r in rows), "rows_total": len(rows)},
        {"obligation": "global Z12 compatibility", "rows_passing": sum(r["obligations"]["global_Z12_compatible"] for r in rows), "rows_total": len(rows)},
        {"obligation": "noncircularity", "rows_passing": sum(r["obligations"]["noncircular_for_target"] for r in rows), "rows_total": len(rows)},
        {"obligation": "no closed-lane replay", "rows_passing": sum(r["obligations"]["no_closed_lane_replay"] for r in rows), "rows_total": len(rows)},
    ]
    theorem = {
        "name": "P3143_T1_vlift_weight_scale_source_obstruction_matrix",
        "statement": "The audited repo-backed scalar candidates provide positive receiver values in 16/20 target rows, but none supplies the full strict-source package for any V_lift coefficient: strict source law, unit-bearing normalization, global Z12 compatibility, noncircularity, and no closed-lane replay.  Therefore P3142's weights/scale remain axiom-branch parameters on current artifacts.",
        "finite_counts": {
            "candidates_tested": len(candidates),
            "targets_tested": len(TARGETS),
            "candidate_target_rows": len(rows),
            "positive_value_rows": sum(r["obligations"]["positive_value"] for r in rows),
            "accepted_strict_source_rows": len(accepted),
        },
    }
    return {
        "status": "P3143_VLIFT_WEIGHT_SCALE_SOURCE_AUDIT_BOUNDED_NO_GO",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "repo_context": {name: inputs[name].get("status") for name in inputs},
        "constructed_object": {
            "name": "Omega_VL source-obligation matrix for V_lift weights/scale",
            "targets": TARGETS,
            "candidate_count": len(candidates),
            "classification": "strict_source_audit_matrix_for_axiom_branch_parameters",
        },
        "candidate_sources": candidates,
        "candidate_target_rows": rows,
        "source_defect_table": source_defect_table,
        "finite_theorem": theorem,
        "decision": {
            "bounded_result": "P3143 performs the narrow source audit requested by P3142.  Existing positive scalar lanes can provide useful receiver magnitudes, but every row fails at least one strict obligation, and no row exports a strict unit-bearing source for mu, w_theta, w_s, or kappa.",
            "negative_export_flags": {
                "vlift_weight_scale_strict_source_exported": False,
                "unit_bearing_selector_action_exported": False,
                "global_Z12_origin_derived": False,
                "strict_QW_2191_discharged": False,
                "strict_selector_closure_exported": False,
                "bridge_completion_exported": False,
                "legacy_role_transfer_exported": False,
                "L_total_exported": False,
                "ToE_closure_exported": False,
            },
            "next_honest_step": "Do not keep recycling scalar magnitude sources for V_lift.  The next honest proof-grade step is exactly one new typed unit-measure object Upsilon_sel that couples the selector local chart to an action measure without importing A_origin/A_lambda; otherwise preserve the P3140-P3143 non-strict axiom-branch/no-strict-source boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    th = payload["finite_theorem"]
    lines = [
        "# P3143/S2093 V_lift weight/scale strict-source audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Constructed object",
        f"- `{payload['constructed_object']['name']}`",
        f"- Targets: `{', '.join(payload['constructed_object']['targets'])}`",
        f"- Classification: `{payload['constructed_object']['classification']}`",
        "",
        "## Finite theorem",
        f"`{th['name']}`: {th['statement']}",
        "",
        "## Finite counts",
    ]
    for key, value in th["finite_counts"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Source-defect table"])
    for row in payload["source_defect_table"]:
        lines.append(f"- `{row['obligation']}`: `{row['rows_passing']}/{row['rows_total']}` rows pass")
    lines.extend(["", "## Candidate summary"])
    for candidate in payload["candidate_sources"]:
        lines.append(f"- `{candidate['id']}` ({candidate['repo_backing']}): value `{candidate['value']}`, strict source `{candidate['strict_source_law']}`, unit-bearing `{candidate['unit_bearing']}`")
    lines.extend(["", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3143/S2093 V_lift weight/scale strict-source audit", "## P3143/S2093 V_lift weight/scale strict-source audit\n\n`P3143/S2093` constructs `Omega_VL`, a candidate-target source-obligation matrix for the `V_lift` parameters `mu`, `w_theta`, `w_s`, and `kappa`.  It audits `5` repo-backed scalar/receiver lanes over `20` target rows: `16/20` rows have positive scalar values, but `0/20` rows export the full strict package of source law, unit-bearing normalization, global `Z12` compatibility, noncircularity, and no closed-lane replay.  Thus `V_lift` weights/scale remain non-strict axiom-branch parameters; no selector closure, bridge completion, role transfer, `L_total`, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3143/S2093 V_lift weights remain unsourced", "## P3143/S2093 V_lift weights remain unsourced\n\n`P3143/S2093` checks whether existing positive scalar lanes can source `V_lift` weights/scale.  The finite obstruction matrix accepts no row as a strict unit-bearing source, so the local variational lift remains a non-strict axiom-branch toy functional rather than a physical action/EOM term.  A next move must introduce a genuinely new unit-measure object `Upsilon_sel` or preserve the no-strict-source boundary.\n")
    append_once(AGENTS, "Current V_lift weight/scale source audit guardrail (P3143/S2093, 2026-07-13)", "## Current V_lift weight/scale source audit guardrail (P3143/S2093, 2026-07-13)\n\n- P3143 audits `5` repo-backed candidate scalar/receiver lanes against the `V_lift` parameters `mu`, `w_theta`, `w_s`, and `kappa`, producing `20` candidate-target rows.\n- Positive magnitudes exist (`16/20` rows), but `0/20` rows export a strict unit-bearing source satisfying source law, global `Z12` compatibility, noncircularity, and no closed-lane replay.\n- Do not recycle entropy, `alpha_geo`, beta/Z_beta, `V_lift` Hessian self-normalization, or `D_HL` receiver obstruction as strict weight/scale source, selector closure, bridge completion, role transfer, `L_total`, or ToE closure.\n- Next honest move: introduce exactly one new typed unit-measure object `Upsilon_sel` coupling the selector local chart to an action measure without importing `A_origin/A_lambda`, or preserve the P3140-P3143 non-strict axiom-branch/no-strict-source boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
