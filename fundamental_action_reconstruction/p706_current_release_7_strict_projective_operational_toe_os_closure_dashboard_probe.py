#!/usr/bin/env python3

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

AS_OF = "2026-03-17"

IN_N705 = (
    GENERATED
    / "n705_current_first_strict_projective_operational_toe_os_support_with_invariant_mass_observable_theorem_summary.json"
)
IN_N701 = GENERATED / "n701_current_strict_projective_operational_toe_closure_theorem_summary.json"
IN_N703 = GENERATED / "n703_current_strict_quadratic_mass_proxy_meaning_definition_theorem_summary.json"
IN_F704 = (
    GENERATED
    / "f704_current_strict_invariant_mass_observable_from_diagonal_local_psi_hessian_eigensystem_export_packet_summary.json"
)
IN_P709 = GENERATED / "p709_current_strict_release_7_os_residual_sign_gauge_irrelevance_audit_probe_summary.json"

# Optional, explicitly non-strict external host-matching probes.
IN_P702 = GENERATED / "p702_current_nonstrict_standard_model_host_matching_from_p696_channel_proxy_probe_summary.json"
IN_P704 = GENERATED / "p704_current_nonstrict_standard_model_host_matching_from_f704_h_psi_eigenvalue_proxy_probe_summary.json"

OUT_JSON = GENERATED / "p706_current_release_7_strict_projective_operational_toe_os_closure_dashboard_probe.json"
OUT_SUMMARY = (
    GENERATED / "p706_current_release_7_strict_projective_operational_toe_os_closure_dashboard_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    strict_prereq = [IN_N705, IN_N701, IN_N703, IN_F704, IN_P709]
    missing = [str(p.relative_to(REPO)) for p in strict_prereq if not p.exists()]
    if missing:
        artifact: dict[str, Any] = {
            "stage": "P706",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    n705 = load_json(IN_N705)
    n701 = load_json(IN_N701)
    n703 = load_json(IN_N703)
    f704 = load_json(IN_F704)
    p709 = load_json(IN_P709)

    n705_tr = (n705.get("theorem_result") or {}) if isinstance(n705, dict) else {}
    n701_tr = (n701.get("theorem_result") or {}) if isinstance(n701, dict) else {}
    n703_tr = (n703.get("theorem_result") or {}) if isinstance(n703, dict) else {}

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any) -> None:
        ok = actual == expected
        checks.append({"id": check_id, "actual": actual, "expected": expected, "pass": ok})
        if not ok:
            blocking.append(check_id)

    # Strict v3 OS support theorem (N705).
    add_check("N705_discharged", bool(n705_tr.get("discharged")), True)
    add_check(
        "N705_os_support_v3_exported",
        bool(n705_tr.get("strict_projective_operational_toe_os_support_packet_v3_exported")),
        True,
    )
    add_check("N705_no_kernel_alone_qw2191_discharge", bool(n705_tr.get("kernel_alone_qw2191_discharge")), False)
    add_check("N705_no_ToE_closure", bool(n705_tr.get("ToE_closure")), False)
    add_check("N705_no_actual_emergent_observer_closure", bool(n705_tr.get("actual_emergent_observer_closure")), False)

    # Operational ToE closure theorem (N701).
    add_check("N701_discharged", bool(n701_tr.get("discharged")), True)
    add_check(
        "N701_operational_toe_closure_strict_projective_os",
        bool(n701_tr.get("operational_toe_closure_strict_projective_os")),
        True,
    )
    add_check("N701_no_kernel_alone_qw2191_discharge", bool(n701_tr.get("kernel_alone_qw2191_discharge")), False)
    add_check("N701_no_ToE_closure", bool(n701_tr.get("ToE_closure")), False)
    add_check(
        "N701_no_actual_emergent_observer_closure",
        bool(n701_tr.get("actual_emergent_observer_closure")),
        False,
    )

    # Quadratic mass-proxy meaning theorem (N703).
    add_check("N703_discharged", bool(n703_tr.get("discharged")), True)
    add_check("N703_quadratic_mass_proxy_meaning_defined", bool(n703_tr.get("quadratic_mass_proxy_meaning_defined")), True)
    add_check("N703_no_physical_unit_identification", bool(n703_tr.get("physical_unit_identification")), False)
    add_check("N703_no_standard_model_host_matching", bool(n703_tr.get("standard_model_host_matching")), False)
    add_check("N703_no_kernel_alone_qw2191_discharge", bool(n703_tr.get("kernel_alone_qw2191_discharge")), False)
    add_check(
        "N703_no_directed_sign_sensitive_physical_orientation_claim",
        bool(n703_tr.get("directed_sign_sensitive_physical_orientation_claim")),
        False,
    )

    # Basis-invariant mass observable export (F704 packet summary).
    add_check("F704_mass_observable_exported", f704.get("status"), "PASS_EXPORTED_STRICT_INVARIANT_MASS_OBSERVABLE_OBJECT")

    # Residual-sign gauge-irrelevance audit (P709).
    add_check(
        "P709_residual_sign_gauge_irrelevance_audited",
        p709.get("status"),
        "PASS_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_AUDITED",
    )

    strict_ok = len(blocking) == 0
    if strict_ok:
        status = "PASS_RELEASE_7_STRICT_PROJECTIVE_OPERATIONAL_OS_CLOSURE_DASHBOARD_READY"
    else:
        status = "P706_REQUIRES_REVIEW_CHANGED_OR_INCOMPLETE_RELEASE_7_STRICT_PROJECTIVE_OPERATIONAL_OS_CLOSURE_STATE"

    nonstrict: dict[str, Any] = {"available": False}
    if IN_P702.exists() or IN_P704.exists():
        nonstrict["available"] = True
        nonstrict["P702"] = load_json(IN_P702) if IN_P702.exists() else None
        nonstrict["P704"] = load_json(IN_P704) if IN_P704.exists() else None

    artifact = {
        "stage": "P706",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "strict_projective_operational_toe_os_closure_dashboard_only",
        "strict_components": {
            "N705": str(IN_N705.relative_to(REPO)),
            "N701": str(IN_N701.relative_to(REPO)),
            "N703": str(IN_N703.relative_to(REPO)),
            "F704": str(IN_F704.relative_to(REPO)),
            "P709": str(IN_P709.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "strict_projective_operational_os_closure_ready": bool(strict_ok),
        "nonstrict_external_host_matching": nonstrict,
        "hard_limits": [
            "no_kernel_alone_global_QW2191_discharge",
            "no_directed_sign_sensitive_physical_orientation_claim",
            "no_standard_model_host_matching_claim_in_strict_scope",
            "no_ToE_closure",
            "no_actual_emergent_observer_closure",
        ],
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
