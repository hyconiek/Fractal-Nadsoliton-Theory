#!/usr/bin/env python3
"""P2029 S979: strict Task-1 renormalization quotient ledger update.

Professor-level decision after P2028:

Do not continue treating Task-1 as "compute four independent coefficients on
scalar B1".  P2028 proves only a rank-3 quotient class modulo the GB/topological
null direction.  This script exports a ledger patch that can be consumed by
later ToE status dashboards without silently promoting the quotient pass into
four-channel uniqueness or tensor-level renormalization closure.
"""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2029_s979_strict_task1_renormalization_quotient_ledger_update.json"
MD = GEN / "p2029_s979_strict_task1_renormalization_quotient_ledger_update.md"

SCHEMA_VERSION = "p2029_s979_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
TASK_ID = 1
TASK_NAME = "Renormalization exact divergence cancellation coefficients"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def find_task(rows: list[dict[str, Any]], task_id: int) -> dict[str, Any]:
    for row in rows:
        if int(row.get("task_id", -1)) == task_id:
            return row
    return {}


def as_bool(value: Any) -> bool:
    return bool(value is True)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2025 = load("p2025_s975_strict_cutkosky_same_scheme_cohomology_amplitude_bridge_seed.json")
    p2027 = load("p2027_s977_strict_b1_rank3_gb_null_adaptive_quadrature_witness.json")
    p2028 = load("p2028_s978_strict_b1_gb_quotient_counterterm_identifiability_theorem.json")

    p2025_gap_task1 = find_task(p2025.get("toe_closure_gaps_7tasks") or [], TASK_ID)
    p2025_numeric_task1 = find_task(p2025.get("task_numeric_evidence_7") or [], TASK_ID)
    p2028_checks = p2028.get("gatekeeper_checks") or {}
    p2027_checks = p2027.get("gatekeeper_checks") or {}
    p2028_q = (p2028.get("quotient_space") or {}).get("target_quotient_coefficients") or {}
    tau_inv = p2028.get("tau_family_invariance") or {}

    quotient_class_pass = (
        p2028.get("local_verdict") == "PASS_QUOTIENT_THEOREM_WITH_TRACE"
        and as_bool(p2028_checks.get("exact_quotient_map_valid"))
        and as_bool(p2028_checks.get("tau_family_invariance_pass"))
        and as_bool(p2028_checks.get("minimum_norm_rep_maps_to_same_quotient"))
    )
    adaptive_scalar_projection_pass = (
        p2027.get("local_verdict") == "PASS_RANK3_QUOTIENT_IDENTIFIABLE_ON_SCALAR_B1_WITH_GB_NULL_TRACE"
        and as_bool(p2027_checks.get("adaptive_transform_gate_pass"))
        and as_bool(p2027_checks.get("full_gram_rank_eq_3"))
        and as_bool(p2027_checks.get("full_gram_nullity_eq_1"))
    )
    independent_gb_identified = as_bool(p2028_checks.get("independent_a_GB_identified"))
    tensor_projection_claimed = as_bool(p2028_checks.get("full_tensor_projection_claimed"))

    evidence_slots = [
        {
            "slot": "p1848_direct_scalar_profiles_consumed",
            "pass": as_bool(p2027_checks.get("p1848_direct_profiles_present")),
            "meaning": "P1950/P2027 consume backend B1 scalar profile exports from P1848.",
        },
        {
            "slot": "adaptive_kdd_quadrature_stable",
            "pass": as_bool(p2027_checks.get("adaptive_transform_gate_pass")),
            "meaning": "Kdd**2-bearing entries are stable under endpoint-transform replay.",
        },
        {
            "slot": "rank3_quotient_class_identified",
            "pass": quotient_class_pass,
            "meaning": "P2028 identifies the scalar B1 counterterm class modulo the GB null direction.",
        },
        {
            "slot": "independent_a_GB_identified",
            "pass": independent_gb_identified,
            "meaning": "This must remain false until a tensor-resolved source separates GB.",
        },
        {
            "slot": "tensor_resolved_projection_exported",
            "pass": tensor_projection_claimed,
            "meaning": "This must remain false until actual tensor component profiles are exported.",
        },
    ]
    positive_support_count = sum(1 for row in evidence_slots[:3] if row["pass"])
    open_hard_blockers = [row["slot"] for row in evidence_slots[3:] if not row["pass"]]

    ledger_patch = {
        "old_p2025_task1_missing": p2025_gap_task1.get("missing", []),
        "replacement_missing": [
            "tensor-resolved B1 curvature-operator projection separating the GB/topological direction",
            "or an explicit theorem licensing Task-1 renormalization strictly in the quotient by span(1,-4,1,-1)",
            "global background transport of the quotient counterterm class beyond scalar B1",
        ],
        "replacement_status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "replacement_result_kind": "PARTIAL_TASK1_RANK3_QUOTIENT_PASS__GB_AND_TENSOR_EXTENSION_OPEN",
        "do_not_use_anymore": [
            "backend-computed a_R2/a_Ric2/a_Riem2/a_GB on strict B1 as if all four were independently identified",
            "minimum-norm four-coefficient vector as a physical a_GB measurement",
        ],
        "safe_task1_statement": (
            "Task-1 has a local strict B1 scalar quotient-class pass: "
            "T(a)=(a_R2+a_GB, a_Ric2-4*a_GB, a_Riem2+a_GB) is identified "
            "with target approximately (1,0,0), while independent a_GB and "
            "tensor-resolved four-channel uniqueness remain open."
        ),
    }

    professor_decision = {
        "decision": "CARRY_QUOTIENT_LEDGER_FORWARD_DO_NOT_CHASE_SCALAR_GB_UNIQUENESS",
        "rationale": [
            "P2028 proves that scalar B1 sees only the quotient class modulo n_GB=(1,-4,1,-1).",
            "Trying to identify a_GB inside the same scalar projection would only choose a section, not add evidence.",
            "The next theoretical bottleneck is tensor-resolved projection or an explicit quotient-form renormalization theorem.",
        ],
        "next_route_preference": "tensor_resolved_projection_source_audit_then_operator_level_projection_attempt",
    }

    local_readiness_score_0_1 = float(positive_support_count / len(evidence_slots))
    status = "OPEN_OBSTRUCTION_WITH_TRACE"
    result_kind = "TASK1_LEDGER_UPDATED_TO_RANK3_QUOTIENT_PASS_TENSOR_EXTENSION_OPEN"
    if not quotient_class_pass or not adaptive_scalar_projection_pass:
        result_kind = "TASK1_LEDGER_UPDATE_OPEN_INPUT_OBSTRUCTION"

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2029",
        "stage_id": "S979",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": status,
        "result_kind": result_kind,
        "route": "strict_only",
        "task_id": TASK_ID,
        "task_name": TASK_NAME,
        "strict_lane_assumptions": [
            "strict_kernel_only",
            "no_legacy_transfer",
            "B1_scalar_quotient_class_not_four_channel_uniqueness",
            "GB_topological_direction_not_independently_identified",
        ],
        "depends_on": {
            "p2025_present": p2025.get("_missing") is None,
            "p2027_present": p2027.get("_missing") is None,
            "p2028_present": p2028.get("_missing") is None,
            "p2028_local_verdict": p2028.get("local_verdict"),
            "p2027_local_verdict": p2027.get("local_verdict"),
        },
        "input_hashes": {
            "p2025_json_sha256": file_sha256(GEN / "p2025_s975_strict_cutkosky_same_scheme_cohomology_amplitude_bridge_seed.json"),
            "p2027_json_sha256": file_sha256(GEN / "p2027_s977_strict_b1_rank3_gb_null_adaptive_quadrature_witness.json"),
            "p2028_json_sha256": file_sha256(GEN / "p2028_s978_strict_b1_gb_quotient_counterterm_identifiability_theorem.json"),
        },
        "professor_decision": professor_decision,
        "p2025_task1_previous_state": {
            "gap_row": p2025_gap_task1,
            "numeric_evidence_row": p2025_numeric_task1,
        },
        "task1_quotient_status": {
            "quotient_class_pass": quotient_class_pass,
            "adaptive_scalar_projection_pass": adaptive_scalar_projection_pass,
            "target_quotient_coefficients": p2028_q,
            "tau_family_max_full_system_residual_linf": tau_inv.get("max_full_system_residual_linf"),
            "tau_family_numerical_tolerance": tau_inv.get("numerical_tolerance"),
            "independent_a_GB_identified": independent_gb_identified,
            "tensor_resolved_projection_exported": tensor_projection_claimed,
            "open_hard_blockers": open_hard_blockers,
        },
        "task1_ledger_patch": ledger_patch,
        "evidence_slots": evidence_slots,
        "local_readiness_score_0_1": local_readiness_score_0_1,
        "honest_verdict": "OPEN_OBSTRUCTION_WITH_TRACE",
        "gatekeeper_checks": {
            "quotient_class_pass": quotient_class_pass,
            "adaptive_scalar_projection_pass": adaptive_scalar_projection_pass,
            "p2025_task1_still_open": p2025_gap_task1.get("status") == "OPEN_OBSTRUCTION_WITH_TRACE",
            "independent_a_GB_not_identified": not independent_gb_identified,
            "tensor_projection_not_claimed": not tensor_projection_claimed,
            "no_toe_closure_claimed": True,
        },
        "theorem_scope": {
            "licensed_statement": ledger_patch["safe_task1_statement"],
            "not_licensed": [
                "four independent scalar B1 counterterm coefficients",
                "a physical measurement of a_GB from the minimum-norm representative",
                "tensor-resolved curvature-operator projection",
                "background-global renormalization",
                "BRST/Cutkosky unitarity closure",
                "QW-2191 selector closure",
                "ToE closure",
            ],
        },
        "false_pass_guard": "This packet updates Task-1 status language only: rank-3 quotient pass is real, but four-channel uniqueness and tensor extension remain open.",
        "next_honest_step": "Build P2030 as a tensor-resolved projection source audit: enumerate what p1848 exports, what component-level tensor data is missing, and the minimal object needed to test GB separation beyond scalar B1.",
        "lay_explanation": "Task-1 nie jest juz blokowany przez brak lokalnych profili B1. Teraz wiemy, ze B1 widzi trzy kombinacje wspolczynnikow, a nie osobny czwarty kierunek GB. Trzeba albo pracowac w tym ilorazie, albo dostarczyc pelna projekcje tensorowa.",
        "environment": {
            "python": platform.python_version(),
        },
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    md = f"""# P2029 S979 Strict Task-1 Renormalization Quotient Ledger Update

Status: `{status}`

Result kind: `{result_kind}`

## Professor Decision

`{professor_decision['decision']}`

P2028 shows that scalar B1 identifies the quotient class

`T(a)=(a_R2+a_GB, a_Ric2-4*a_GB, a_Riem2+a_GB)`

modulo `n_GB=(1,-4,1,-1)`.  Therefore the old Task-1 wording from P2025,
`backend-computed a_R2/a_Ric2/a_Riem2/a_GB on strict B1`, is superseded.

## Updated Task-1 State

Local quotient class pass: `{quotient_class_pass}`.

Adaptive scalar projection pass: `{adaptive_scalar_projection_pass}`.

Target quotient coefficients:
- R2_bar={float(p2028_q.get('R2_bar', 0.0)):.16e}
- Ric2_bar={float(p2028_q.get('Ric2_bar', 0.0)):.16e}
- Riem2_bar={float(p2028_q.get('Riem2_bar', 0.0)):.16e}

Hard blockers that remain:
- independent `a_GB` is not identified;
- tensor-resolved projection is not exported;
- background-global transport of the quotient class is not yet proven.

## Replacement Missing List

1. tensor-resolved B1 curvature-operator projection separating the GB/topological direction;
2. or an explicit theorem licensing Task-1 renormalization strictly in the quotient by `span(1,-4,1,-1)`;
3. global background transport of the quotient counterterm class beyond scalar B1.

## Honest Interpretation

Task-1 has progressed from "no direct backend profiles" to a real local
rank-3 quotient theorem.  It has not become four-channel renormalization
closure.  The next research move should be a tensor-resolved source audit and
operator-level projection attempt, not another scalar GB coefficient fit.
"""
    MD.write_text(md, encoding="utf-8")
    print(OUT)
    print(MD)


if __name__ == "__main__":
    main()
