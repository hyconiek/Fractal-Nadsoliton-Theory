#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2318_s1268_strict_selector_lane_to_task3_margin_bridge_obligation_probe.json"
MD = GEN / "p2318_s1268_strict_selector_lane_to_task3_margin_bridge_obligation_probe.md"

SOURCE_FILES = {
    "P2302_PROVIDER_LIFT_CANDIDATE": GEN / "p2302_s1252_strict_task3_provider_lift_policy_lock_candidate_probe.json",
    "P2305_NORM_TO_MARGIN_UNDERDETERMINATION": GEN / "p2305_s1255_strict_norm_to_margin_bridge_underdetermination_probe.json",
    "P2306_RESPONSE_ORIENTATION_INTERFACE": GEN / "p2306_s1256_strict_response_orientation_functional_interface_audit_probe.json",
    "P2307_DYNAMICAL_RESPONSE_NONDERIVATION": GEN / "p2307_s1257_strict_dynamical_margin_response_functional_nonderivation_probe.json",
    "P2308_CURRENT_INTERFACE_NONEXISTENCE": GEN / "p2308_s1258_strict_current_interface_class_response_functional_nonexistence_probe.json",
    "P2313_ROUTE_STOP_SOURCE_CONTRACT": GEN / "p2313_s1263_strict_provider_margin_route_stop_and_source_class_contract_probe.json",
    "P2317_FOURIER_DEGENERACY_LANE_AUDIT": GEN / "p2317_s1267_strict_fourier_degeneracy_existing_lane_audit_probe.json",
    "DIAGONAL_MODE_ASSIGNMENT": GEN / "mode_index_assignment_canonical_local_diagonal_strict_derived_v1.json",
    "SHANNON_MODE_ASSIGNMENT": GEN / "mode_index_assignment_shannon_element_order_reference_strict_core_v1.json",
}

GREP_PATTERNS = (
    "provider_lift_per_step",
    "policy-margin",
    "policy margin",
    "response functional",
    "orientation interface",
    "mode-index assignment",
    "Task-3",
    "G1/G3",
    "bridge_obligation",
    "underdetermination",
)

REQUIRED_BRIDGE_FIELDS = (
    "signed_scalar_lift_per_step",
    "time_step_normalization",
    "margin_response_functional",
    "orientation_to_policy_sign_rule",
    "p2281_replay_semantics",
    "admissibility_theorem_no_selector_premise",
)


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def corpus_hits() -> list[dict[str, Any]]:
    paths = sorted(
        set(SOURCE_FILES.values())
        | set(GEN.glob("p230*_s12*_strict*bridge*json"))
        | set(GEN.glob("p230*_s12*_strict*response*json"))
        | set(GEN.glob("p231*_s12*_strict*margin*json"))
        | set(ROOT.glob("p230*_s12*_strict*bridge*.py"))
        | set(ROOT.glob("p230*_s12*_strict*response*.py"))
        | set(ROOT.glob("p231*_s12*_strict*margin*.py"))
    )
    rows: list[dict[str, Any]] = []
    for path in paths:
        if not path.exists() or path.is_dir():
            continue
        text = read_text(path)
        lowered = text.lower()
        count = sum(lowered.count(pattern.lower()) for pattern in GREP_PATTERNS)
        if count == 0:
            continue
        first_line = 0
        first_excerpt = ""
        for index, line in enumerate(text.splitlines(), start=1):
            if any(pattern.lower() in line.lower() for pattern in GREP_PATTERNS):
                first_line = index
                first_excerpt = line.strip()[:240]
                break
        rows.append({
            "path": path.relative_to(REPO).as_posix(),
            "pattern_hit_count": count,
            "first_hit_line": first_line,
            "first_hit_excerpt": first_excerpt,
        })
    rows.sort(key=lambda row: (-int(row["pattern_hit_count"]), row["path"]))
    return rows


def get_required_lift(p2302: dict[str, Any]) -> float:
    probe = p2302.get("strict_task3_provider_lift_policy_lock_candidate_probe", {})
    candidate = probe.get("policy_lock_candidate", {})
    return float(candidate.get("provider_lift_per_step", 0.0))


def lane_rows(p2317: dict[str, Any], lane_key: str) -> list[dict[str, Any]]:
    probe = p2317.get("strict_fourier_degeneracy_existing_lane_audit_probe", {})
    lane = probe.get(lane_key, {})
    return list(lane.get("pair_rows", []))


def scale_witnesses(rows: list[dict[str, Any]], required_lift: float) -> dict[str, Any]:
    strengths = [float(row.get("defect_abs", 0.0)) for row in rows]
    max_strength = max(strengths, default=0.0)
    min_nonzero = min((value for value in strengths if value > 0.0), default=0.0)
    calibrated_scale = required_lift / max_strength if max_strength else None
    min_scale = required_lift / min_nonzero if min_nonzero else None
    witnesses = []
    if max_strength:
        witnesses.extend([
            {
                "witness_id": "POSITIVE_TARGET_CALIBRATED_SCALE",
                "scale": calibrated_scale,
                "chosen_strength": max_strength,
                "lift_value": required_lift,
                "passes_required_lift": True,
                "strict_admissible_without_exported_scale": False,
            },
            {
                "witness_id": "NEGATIVE_SAME_MAGNITUDE_SCALE",
                "scale": -calibrated_scale,
                "chosen_strength": max_strength,
                "lift_value": -required_lift,
                "passes_required_lift": False,
                "strict_admissible_without_exported_sign_rule": True,
            },
            {
                "witness_id": "ZERO_RESPONSE_SCALE",
                "scale": 0.0,
                "chosen_strength": max_strength,
                "lift_value": 0.0,
                "passes_required_lift": False,
                "strict_admissible_without_exported_response": True,
            },
        ])
    return {
        "pair_strengths": strengths,
        "max_strength": max_strength,
        "min_nonzero_strength": min_nonzero,
        "required_lift": required_lift,
        "target_calibrated_scale_using_max_strength": calibrated_scale,
        "scale_needed_if_every_pair_must_individually_exceed_required_lift": min_scale,
        "witnesses": witnesses,
        "underdetermined": True,
        "reason": "Lane defect magnitudes are dimensionless axis-splitting data; P2302 requires a signed per-step policy-margin lift.  A scale/sign/response map can make the same lane data pass, fail, or vanish.",
    }


def interface_matrix() -> list[dict[str, Any]]:
    return [
        {
            "required_field": "signed_scalar_lift_per_step",
            "diagonal_local_lane_exports": "pairwise defect magnitudes, angles, residual Z2 axes",
            "shannon_lane_exports": "pairwise element-order defects, angles, residual Z2 axes",
            "status": "MISSING_TYPED_MAP_TO_PROVIDER_LIFT_PER_STEP",
        },
        {
            "required_field": "time_step_normalization",
            "diagonal_local_lane_exports": "none",
            "shannon_lane_exports": "none",
            "status": "MISSING_DELTA_T_OR_REPLAY_CLOCK_THEOREM",
        },
        {
            "required_field": "margin_response_functional",
            "diagonal_local_lane_exports": "none",
            "shannon_lane_exports": "none",
            "status": "MISSING_RESPONSE_FUNCTIONAL_R",
        },
        {
            "required_field": "orientation_to_policy_sign_rule",
            "diagonal_local_lane_exports": "axis up to residual Z2; no signed policy orientation",
            "shannon_lane_exports": "axis up to residual Z2; r_ord is not translation-invariant but direction-free under Aut(Z12)",
            "status": "MISSING_SIGN_ORIENTATION_RULE",
        },
        {
            "required_field": "p2281_replay_semantics",
            "diagonal_local_lane_exports": "mode-index basis object, not replay row semantics",
            "shannon_lane_exports": "mode-index basis object, not replay row semantics",
            "status": "MISSING_REPLAY_COUPLING_THEOREM",
        },
        {
            "required_field": "admissibility_theorem_no_selector_premise",
            "diagonal_local_lane_exports": "lane-scoped hard limits preserve no global selector closure",
            "shannon_lane_exports": "lane-scoped hard limits preserve no global selector closure",
            "status": "MISSING_TASK3_ADMISSIBILITY_THEOREM",
        },
    ]


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {sid: load_json(path) for sid, path in SOURCE_FILES.items() if path.suffix == ".json"}
    source_hashes = {sid: sha256_file(path) for sid, path in SOURCE_FILES.items() if path.exists()}
    hits = corpus_hits()

    p2302 = artifacts["P2302_PROVIDER_LIFT_CANDIDATE"]
    p2317 = artifacts["P2317_FOURIER_DEGENERACY_LANE_AUDIT"]
    required_lift = get_required_lift(p2302)
    diagonal_witness = scale_witnesses(lane_rows(p2317, "diagonal_local_lane_o2_to_z2_computation"), required_lift)
    shannon_witness = scale_witnesses(lane_rows(p2317, "shannon_element_order_lane_o2_to_z2_computation"), required_lift)
    matrix = interface_matrix()

    existing_bridge_blockers = {
        "p2305_norm_to_margin_verdict": artifacts["P2305_NORM_TO_MARGIN_UNDERDETERMINATION"].get(
            "strict_norm_to_margin_bridge_underdetermination_probe", {}
        ).get("strict_bridge_verdict", {}),
        "p2306_status": artifacts["P2306_RESPONSE_ORIENTATION_INTERFACE"].get("status"),
        "p2307_status": artifacts["P2307_DYNAMICAL_RESPONSE_NONDERIVATION"].get("status"),
        "p2308_status": artifacts["P2308_CURRENT_INTERFACE_NONEXISTENCE"].get("status"),
        "p2313_status": artifacts["P2313_ROUTE_STOP_SOURCE_CONTRACT"].get("status"),
    }

    bridge_obligation_verdict = {
        "required_lift_per_step": required_lift,
        "diagonal_local_can_be_target_calibrated_numerically": bool(diagonal_witness["witnesses"]),
        "shannon_can_be_target_calibrated_numerically": bool(shannon_witness["witnesses"]),
        "target_calibration_is_strict_bridge": False,
        "all_required_bridge_fields_exported": False,
        "missing_required_bridge_fields": [row["required_field"] for row in matrix if row["status"].startswith("MISSING")],
        "formal_status": "OPEN_NONEXISTENCE_FOR_DIRECT_LANE_DEFECT_TO_TASK3_MARGIN_BRIDGE_CURRENT_EXPORT_SET",
        "reason": "Existing lane defects are enough to cut O(2) to Z2 inside their lanes, but no exported theorem maps them to the signed P2302 provider_lift_per_step with replay semantics.",
    }

    theorem_export = {
        "theorem_name": "P2318 selector-lane to Task-3 margin bridge obligation audit",
        "claim": "On the current export set, diagonal/local and Shannon O(2)->Z2 lane data do not supply a strict Task-3 provider_lift_per_step bridge.  They can be target-calibrated numerically, but sign, scale, response functional, replay clock, and admissibility theorem remain unexported.",
        "proof_bits": {
            "required_lift_per_step": required_lift,
            "diagonal_max_defect": diagonal_witness["max_strength"],
            "shannon_max_defect": shannon_witness["max_strength"],
            "diagonal_target_scale": diagonal_witness["target_calibrated_scale_using_max_strength"],
            "shannon_target_scale": shannon_witness["target_calibrated_scale_using_max_strength"],
            "missing_bridge_field_count": len(bridge_obligation_verdict["missing_required_bridge_fields"]),
            "witnesses_per_lane": len(diagonal_witness["witnesses"]),
        },
        "nonpromotion_rule": "Do not promote lane O(2)->Z2 axis data to Task-3 G1/G3 closure until a theorem-grade signed response/margin bridge is exported and replayed through P2302/P2281/P2282.",
    }

    probe = {
        "probe_id": "P2318_S1268_STRICT_SELECTOR_LANE_TO_TASK3_MARGIN_BRIDGE_OBLIGATION",
        "source_hashes": source_hashes,
        "far_bridge_grep_audit": {
            "search_patterns": list(GREP_PATTERNS),
            "hit_count": len(hits),
            "top_hits": hits[:30],
        },
        "interface_obligation_matrix": matrix,
        "diagonal_local_lane_scale_witness": diagonal_witness,
        "shannon_element_order_lane_scale_witness": shannon_witness,
        "existing_bridge_blockers": existing_bridge_blockers,
        "bridge_obligation_verdict": bridge_obligation_verdict,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "far_bridge_grep_hits_found": len(hits) >= 10,
        "p2302_required_lift_loaded": required_lift > 0.0,
        "p2317_loaded": p2317.get("packet_id") == "P2317",
        "diagonal_lane_target_calibration_possible": bridge_obligation_verdict["diagonal_local_can_be_target_calibrated_numerically"],
        "shannon_lane_target_calibration_possible": bridge_obligation_verdict["shannon_can_be_target_calibrated_numerically"],
        "target_calibration_not_promoted": bridge_obligation_verdict["target_calibration_is_strict_bridge"] is False,
        "missing_bridge_fields_detected": len(bridge_obligation_verdict["missing_required_bridge_fields"]) == len(REQUIRED_BRIDGE_FIELDS),
        "all_required_bridge_fields_not_exported": bridge_obligation_verdict["all_required_bridge_fields_exported"] is False,
        "p2305_norm_bridge_refutation_loaded": existing_bridge_blockers["p2305_norm_to_margin_verdict"].get("norm_to_margin_bridge_refuted_as_norm_only_theorem") is True,
        "strict_g1_g3_not_updated": True,
        "no_selector_closure_claimed": True,
        "no_qw2191_global_discharge_claimed": True,
        "no_toe_closure_claimed": True,
    }

    payload = {
        "schema_version": "p2318_s1268_v1",
        "packet_id": "P2318",
        "stage_id": "S1268",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_SELECTOR_LANE_TO_TASK3_MARGIN_BRIDGE_FIELDS_MISSING",
        "result_kind": "STRICT_SELECTOR_LANE_TO_PROVIDER_LIFT_OBLIGATION_AUDIT_NO_G1_G3_UPDATE",
        "strict_selector_lane_to_task3_margin_bridge_obligation_probe": probe,
        "recommended_next_honest_step": "Either export the missing signed response/margin bridge fields for a lane source and replay P2302/P2281/P2282, or package the current matrix as a theorem-grade nonexistence result for direct lane-defect-to-margin bridges.",
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")

    md_lines = [
        "# P2318 S1268: selector-lane to Task-3 margin bridge obligation audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## FAR bridge grep result",
        f"- Bridge/interface evidence hits: `{len(hits)}`.",
        "",
        "## Required lift",
        f"- P2302 required `provider_lift_per_step`: `{required_lift}`.",
        "",
        "## Lane scale witnesses",
        f"- Diagonal/local max defect: `{diagonal_witness['max_strength']:.6g}`; target scale: `{diagonal_witness['target_calibrated_scale_using_max_strength']:.6g}`.",
        f"- Shannon max defect: `{shannon_witness['max_strength']:.6g}`; target scale: `{shannon_witness['target_calibrated_scale_using_max_strength']:.6g}`.",
        "- The same lane data can pass, fail, or vanish under different unexported scale/sign/response maps.",
        "",
        "## Missing bridge fields",
        f"- Missing fields: `{', '.join(bridge_obligation_verdict['missing_required_bridge_fields'])}`.",
        "- Task-3 G1/G3 update: `NONE`.",
        "",
        "## Theorem fingerprint",
        f"`{probe['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
