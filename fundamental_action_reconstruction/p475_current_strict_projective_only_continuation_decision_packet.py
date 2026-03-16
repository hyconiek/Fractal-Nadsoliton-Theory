#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

PROJECTIVE_STATE = GENERATED / "selector_state_global_c_v1_projective_strict_v1.json"
P474_GLUING_SUMMARY = (
    GENERATED
    / "p474_current_strict_global_projective_selector_state_gluing_consistency_audit_probe_summary.json"
)
H37_SIGN_AUDIT_SUMMARY = GENERATED / "h37_sign_distinction_state_audit_summary.json"

N502_THEOREM = (
    ROOT
    / "N502_CURRENT_FIRST_STRICT_RESIDUAL_Z2_SIGN_FREEZE_FOR_EXPORTED_DOWNSTREAM_OBJECTS_GAUGE_IRRELEVANCE_PACKAGE_THEOREM.md"
)
N519_THEOREM = (
    ROOT
    / "N519_CURRENT_FIRST_STRICT_RESIDUAL_Z2_SIGN_FREEZE_EXTENDED_TO_GLOBAL_PROJECTIVE_SELECTOR_OBJECTS_GAUGE_IRRELEVANCE_THEOREM.md"
)
T171_TARGET = ROOT / "T171_CURRENT_STRICT_GLOBAL_DIRECTED_SELECTOR_STATE_DATUM_TARGET_SPEC.md"

P11_FACTORIZATION_SUMMARY = (
    GENERATED
    / "p11_existing_kernel_feedback_to_explicit_h3_chain_factorization_map_probe_summary.json"
)

OUT = (
    GENERATED / "p475_current_strict_projective_only_continuation_decision_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p475_current_strict_projective_only_continuation_decision_packet_summary.json"
)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _p(path: Path) -> str:
    try:
        return str(path.relative_to(REPO))
    except ValueError:
        return str(path)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = {
        "projective_state": PROJECTIVE_STATE,
        "p474_gluing_summary": P474_GLUING_SUMMARY,
        "h37_sign_audit_summary": H37_SIGN_AUDIT_SUMMARY,
        "N502_theorem": N502_THEOREM,
        "N519_theorem": N519_THEOREM,
        "T171_target_spec": T171_TARGET,
    }
    missing_required = [k for k, p in required.items() if not p.is_file()]
    if missing_required:
        payload = {
            "stage": "P475",
            "date": "2026-03-16",
            "goal": "declare_professorial_projective_only_continuation_for_strict_global_selector_state__defer_directed_sign_lift",
            "status": "FAIL_MISSING_REQUIRED_INPUTS",
            "missing_required_inputs": missing_required,
            "required_paths": {k: _p(v) for k, v in required.items()},
            "no_false_pass": True,
            "hard_limits": [
                "no theorem-level pass",
                "no full-closure pass",
                "no claim that H37/T171 is discharged",
                "no sign-sensitive/directed selector state claim",
                "no claim of global QW-2191 discharge",
                "no claim of ToE closure",
            ],
        }
        OUT.write_text(
            json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
        )
        OUT_SUMMARY.write_text(
            json.dumps(
                {
                    "stage": payload["stage"],
                    "status": payload["status"],
                    "decision": None,
                    "recommended_next_strict_target": "H37",
                    "no_false_pass": True,
                },
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    projective = _read_json(PROJECTIVE_STATE)
    p474 = _read_json(P474_GLUING_SUMMARY)
    h37 = _read_json(H37_SIGN_AUDIT_SUMMARY)

    sign_gauge = (
        (projective.get("state_type") or {}).get("sign_gauge")
        if isinstance(projective.get("state_type"), dict)
        else None
    )
    p474_pass = bool(p474.get("overall_pass"))

    if projective.get("no_false_pass") is not True:
        raise ValueError(
            "expected no_false_pass=true in the global projective selector state object"
        )
    if p474.get("no_false_pass") is not True:
        raise ValueError("expected no_false_pass=true in P474 summary")
    if not sign_gauge:
        raise ValueError(
            "expected state_type.sign_gauge=true in the global projective selector state object"
        )
    if not p474_pass:
        raise ValueError("expected overall_pass=true in P474 summary")

    # P475 is a *decision*; we still validate the current strict status boundary to avoid a false pass.
    strict_sign_sensitive_present = False
    if isinstance(h37, dict):
        strict_sign_sensitive_present = (
            "SIGN_SENSITIVE" in str(h37.get("status") or "").upper()
        )

    p11_next: dict[str, Any] | None = None
    if P11_FACTORIZATION_SUMMARY.exists():
        try:
            p11 = _read_json(P11_FACTORIZATION_SUMMARY)
            p11_next = {
                "stage": "P11",
                "status": p11.get("status"),
                "required_next_step": p11.get("required_next_step"),
                "missing_upstream_objects": p11.get("missing_upstream_objects"),
                "summary_path": _p(P11_FACTORIZATION_SUMMARY),
            }
        except Exception:
            p11_next = {
                "stage": "P11",
                "summary_path": _p(P11_FACTORIZATION_SUMMARY),
                "parse_error": True,
            }

    status = "PASS_DECISION_DECLARED_PROJECTIVE_ONLY_CONTINUATION_SELECTED"
    decision = "PROJECTIVE_ONLY_CONTINUATION_SELECTED"

    payload = {
        "stage": "P475",
        "date": "2026-03-16",
        "goal": "declare_professorial_projective_only_continuation_for_strict_global_selector_state__defer_directed_sign_lift",
        "status": status,
        "decision": decision,
        "decision_basis": {
            "projective_state_object_present": _p(PROJECTIVE_STATE),
            "projective_state_sign_gauge": bool(sign_gauge),
            "p474_gluing_audit_pass": bool(p474_pass),
            "p474_summary": _p(P474_GLUING_SUMMARY),
            "gauge_irrelevance_packages": {
                "N502": _p(N502_THEOREM),
                "N519": _p(N519_THEOREM),
            },
            "directed_frontier_target_spec": _p(T171_TARGET),
            "strict_sign_sensitive_datum_present": bool(strict_sign_sensitive_present),
            "strict_boundary_note": "H37/T171 remain open for directed continuation; no strict sign-sensitive observable is promoted by this decision.",
        },
        "continuation": {
            "selected": "projective_only",
            "meaning": "treat exported global selector state as a projective/ray-level physical state object for the declared strict closure stack; keep residual sign as gauge/convention where proven irrelevant; defer directed/sign-sensitive lift attempts to an explicit future branch.",
            "directed_branch_status": "OPEN",
        },
        "recommended_next_strict_target": {
            "target": "P11",
            "note": "After accepting projective-only continuation, proceed with strict-only ToE closure tasks that depend only on projectors/spans. The current missing object on the existing-kernel-feedback → K_obs route remains the explicit factorization/equivalence map into the explicit H3 chain (P10/P11).",
            "p11_next": p11_next,
        },
        "hard_limits": [
            "no theorem-level pass",
            "no full-closure pass",
            "no claim that H37/T171 is discharged",
            "no sign-sensitive/directed selector state claim",
            "no claim of global QW-2191 discharge",
            "no claim of ToE closure",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": payload["stage"],
        "status": payload["status"],
        "decision": payload["decision"],
        "selected_continuation": payload["continuation"]["selected"],
        "recommended_next_strict_target": payload["recommended_next_strict_target"][
            "target"
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    print(OUT)


if __name__ == "__main__":
    main()

