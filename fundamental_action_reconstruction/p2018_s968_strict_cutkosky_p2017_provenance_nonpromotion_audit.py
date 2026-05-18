#!/usr/bin/env python3
"""P2018 S968 strict Cutkosky P2017 provenance non-promotion audit.

Professorial corrective step after P2017: do not add another amplitude ansatz.
Instead, audit the P2017 quadrature-candidate witness against the earlier P1953
dressed-amplitude interface contract and mechanically enforce that it cannot be
promoted to backend-derived Cutkosky closure while contract fields are missing.
"""
from __future__ import annotations

import json
import platform
from pathlib import Path
from typing import Any

import numpy as np
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2018_s968_strict_cutkosky_p2017_provenance_nonpromotion_audit.json"
TS = "2026-05-18T00:00:00+00:00"
P2017_NAME = "p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.json"
P1953_NAME = "p1953_s903_strict_dressed_cutkosky_amplitude_availability_audit.json"

REQUIRED_P2017_OPEN_RESULT = "OPEN_PROVENANCE_GAP_WITH_STRICT_QUADRATURE_TRACE"
REQUIRED_P2017_STATUS = "OPEN_OBSTRUCTION_WITH_TRACE"


def load(name: str) -> dict[str, Any]:
    path = GEN / name
    if not path.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def bool_false_keys(mapping: dict[str, Any]) -> list[str]:
    return [key for key, value in mapping.items() if value is False]


def has_key_deep(obj: Any, wanted: str) -> bool:
    if isinstance(obj, dict):
        return wanted in obj or any(has_key_deep(value, wanted) for value in obj.values())
    if isinstance(obj, list):
        return any(has_key_deep(value, wanted) for value in obj)
    return False


def contract_row(field: str, verdict: str, evidence: str) -> dict[str, str]:
    return {"required_field": field, "verdict": verdict, "evidence": evidence}


def evaluate_p1953_contract(p2017: dict[str, Any], required_fields: list[str]) -> list[dict[str, str]]:
    """Evaluate P2017 against the P1953 dressed-amplitude interface contract.

    This deliberately allows partial evidence only where P2017 actually exports a
    compatible object.  Partial evidence is not enough for promotion.
    """
    rows: list[dict[str, str]] = []
    for field in required_fields:
        if field == "channel = graviton->gauge_gauge":
            ok = p2017.get("channel") == "graviton->gauge_gauge"
            rows.append(contract_row(field, "PASS_PRESENT" if ok else "FAIL_MISSING", "P2017 channel field"))
        elif field == "phase_space_measure and integration_domain":
            contract = p2017.get("quadrature_amplitude_candidate_contract", {})
            has_interval = "phase_space_interval" in contract
            rows.append(contract_row(
                field,
                "PARTIAL_INTERVAL_ONLY" if has_interval else "FAIL_MISSING",
                "P2017 exports an ansatz interval but not a strict-derived phase-space measure.",
            ))
        elif field == "uncertainty_or_exactness_certificate":
            has_quad_error = has_key_deep(p2017, "quad_error_3x3")
            rows.append(contract_row(
                field,
                "PARTIAL_QUADRATURE_ERROR_ONLY" if has_quad_error else "FAIL_MISSING",
                "P2017 exports quadrature errors for candidate tensors, not theorem-grade exactness/uncertainty for a dressed amplitude.",
            ))
        elif field in {
            "M_dressed_common_basis",
            "AbsM_dressed_squared_common_basis",
            "DiscM_common_basis",
            "CutSum_common_basis",
            "DiscM_minus_CutSum_simplified",
            "dressed_graviton_propagator_pole_list",
            "residue_values_per_pole",
            "ghost_sector_exclusion_trace",
        }:
            rows.append(contract_row(field, "FAIL_MISSING", "P2017 is a quadrature-candidate witness and exports no dressed common-basis amplitude object."))
        elif field == "external_state_projectors including BRST physical-state projector":
            rows.append(contract_row(field, "FAIL_MISSING", "No BRST physical-state projector is exported by P2017."))
        elif field == "gauge_fixing_family and xi value":
            rows.append(contract_row(field, "FAIL_MISSING", "P2017 does not declare a gauge-fixing family or xi value."))
        elif field == "scheme = MSbar_B1_seed":
            rows.append(contract_row(field, "FAIL_MISSING", "P2017 does not lock its candidate tensor to MSbar_B1_seed."))
        else:
            rows.append(contract_row(field, "FAIL_UNCLASSIFIED", "No P2018 matcher exists for this required field; promotion is blocked."))
    return rows


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2017 = load(P2017_NAME)
    p1953 = load(P1953_NAME)

    provenance = p2017.get("provenance_gatekeeper_checks", {})
    numerical = p2017.get("numerical_gatekeeper_checks", {})
    combined = p2017.get("gatekeeper_checks", {})
    contract = p1953.get("dressed_amplitude_interface_contract", {})
    required_fields = list(contract.get("required_fields", []))
    contract_rows = evaluate_p1953_contract(p2017, required_fields)
    missing_contract_fields = [
        row["required_field"] for row in contract_rows if row["verdict"].startswith("FAIL")
    ]
    partial_contract_fields = [
        row["required_field"] for row in contract_rows if row["verdict"].startswith("PARTIAL")
    ]

    cov = np.array(
        p2017.get("quadrature_channel_covariance_candidate", {}).get(
            "C4_channel_covariance_from_quadrature_traces", []
        ),
        dtype=float,
    )

    missing_provenance = bool_false_keys(provenance)
    numerical_false = bool_false_keys(numerical)
    provenance_gap_count = len(missing_provenance)

    if cov.shape == (4, 4):
        eigvals = np.linalg.eigvalsh(cov)
        fro_norm = float(la.norm(cov, "fro"))
        numerical_rank = int(np.linalg.matrix_rank(cov, tol=1e-18))
        min_eig = float(np.min(eigvals))
        max_eig = float(np.max(eigvals))
    else:
        eigvals = np.array([], dtype=float)
        fro_norm = float("nan")
        numerical_rank = 0
        min_eig = float("nan")
        max_eig = float("nan")

    p2017_result_open = p2017.get("result_kind") == REQUIRED_P2017_OPEN_RESULT
    p2017_status_open = p2017.get("status") == REQUIRED_P2017_STATUS
    p2017_declares_provenance_gap = combined.get("provenance_gap_declared") is True
    p2017_blocks_false_pass = combined.get("false_pass_blocked_by_provenance_gate") is True
    has_provenance_gap = provenance_gap_count > 0
    contract_complete = bool(required_fields) and not missing_contract_fields and not partial_contract_fields

    nonpromotion_gate = {
        "p1953_contract_present": bool(required_fields),
        "p2017_present": p2017.get("packet_id") == "P2017",
        "p2017_result_kind_open": p2017_result_open,
        "p2017_status_open": p2017_status_open,
        "p2017_declares_provenance_gap": p2017_declares_provenance_gap,
        "p2017_blocks_false_pass": p2017_blocks_false_pass,
        "p2017_has_missing_provenance_gates": has_provenance_gap,
        "p2017_fails_p1953_dressed_amplitude_contract": not contract_complete,
        "numerical_diagnostics_do_not_override_provenance": has_provenance_gap
        and p2017.get("diagnostic_result_kind") == "PASS_STRICT_KERNEL_QUADRATURE_NUMERICS",
        "no_closed_cutkosky_result_exported": p2017.get("result_kind") != "PASS_BACKEND_LOOP_AMPLITUDE_TENSOR_QUADRATURE_WITNESS",
        "covariance_candidate_shape_audited": cov.shape == (4, 4),
        "covariance_candidate_psd_or_missing_nonfatal": bool(cov.shape == (4, 4) and min_eig >= -1e-15),
    }

    promotion_verdict = (
        "BLOCK_PROMOTION_TO_BACKEND_CUTKOSKY_CLOSURE"
        if all(nonpromotion_gate.values())
        else "AUDIT_REQUIRES_MANUAL_REVIEW"
    )

    promote = sp.Symbol("promote_to_backend_cutkosky_closure")
    provenance_complete = sp.Symbol("provenance_complete")
    contract_complete_sym = sp.Symbol("p1953_contract_complete")
    symbolic_rule = sp.Implies(promote, provenance_complete & contract_complete_sym)
    current_rule_value = bool(symbolic_rule.subs({
        promote: True,
        provenance_complete: False,
        contract_complete_sym: contract_complete,
    }))

    out = {
        "ledger_id": "P2018_S968_STRICT_CUTKOSKY_P2017_PROVENANCE_NONPROMOTION_AUDIT",
        "packet_id": "P2018",
        "stage_id": "S968",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "depends_on": {"p2017": p2017.get("ledger_id"), "p1953": p1953.get("packet_id")},
        "audit_target": P2017_NAME,
        "contract_source": P1953_NAME,
        "audit_scope": "nonpromotion audit of P2017 quadrature candidate against the P1953 dressed-amplitude interface; no new amplitude ansatz introduced",
        "p2017_result_summary": {
            "diagnostic_result_kind": p2017.get("diagnostic_result_kind"),
            "result_kind": p2017.get("result_kind"),
            "status": p2017.get("status"),
            "false_pass_guard": p2017.get("false_pass_guard"),
        },
        "p1953_contract_audit": {
            "object_name": contract.get("object_name"),
            "acceptance_rule": contract.get("acceptance_rule"),
            "rows": contract_rows,
            "missing_contract_fields": missing_contract_fields,
            "partial_contract_fields": partial_contract_fields,
            "contract_complete": contract_complete,
        },
        "missing_provenance_gates": missing_provenance,
        "numerical_gate_failures": numerical_false,
        "provenance_gap_count": provenance_gap_count,
        "quadrature_covariance_candidate_audit": {
            "shape": list(cov.shape),
            "frobenius_norm": fro_norm,
            "numerical_rank_tol_1e_minus_18": numerical_rank,
            "eigvals": eigvals.tolist(),
            "min_eig": min_eig,
            "max_eig": max_eig,
        },
        "symbolic_nonpromotion_rule": {
            "rule": "promote_to_backend_cutkosky_closure => (provenance_complete and p1953_contract_complete)",
            "sympy_expression": str(symbolic_rule),
            "current_provenance_complete": False,
            "current_p1953_contract_complete": contract_complete,
            "current_rule_value_if_promotion_attempted": current_rule_value,
            "current_promote_allowed": False,
        },
        "gatekeeper_checks": nonpromotion_gate,
        "promotion_verdict": promotion_verdict,
        "result_kind": "PASS_PROVENANCE_NONPROMOTION_AUDIT" if all(nonpromotion_gate.values()) else "OPEN_AUDIT_REVIEW_REQUIRED",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "P2018 is a non-promotion audit against the P1953 dressed-amplitude contract. It prevents P2017 numerical quadrature diagnostics from being read as backend-derived Cutkosky closure, all-state unitarity, QW-2191 discharge, or ToE closure.",
        "next_honest_step": "Build one strict-side dressed-amplitude interface component from the P1953 contract, starting with a same-scheme BRST-projected M_dressed_common_basis component, then rerun P2018.",
        "toe_progress": "Improves ToE rigor by tying the P2017 non-promotion boundary to the earlier P1953 dressed-amplitude contract instead of relying only on hard-open local provenance flags.",
        "lay_explanation": "To jest kontrola bezpieczeństwa oparta na wcześniejszej liście wymagań: nawet dobre całki P2017 nie spełniają kontraktu prawdziwej ubranej amplitudy, więc nie wolno nazwać ich pełnym dowodem unitarności.",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": __import__("scipy").__version__,
            "sympy": sp.__version__,
        },
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2018] wrote audit: {OUT}")


if __name__ == "__main__":
    main()
