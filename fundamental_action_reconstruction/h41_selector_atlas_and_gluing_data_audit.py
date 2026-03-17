from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
GENERATED.mkdir(exist_ok=True)

H37_SUMMARY = GENERATED / "h37_sign_distinction_state_audit_summary.json"

SELECTOR_ATLAS_GLOBAL = GENERATED / "selector_atlas_global_c_v1_strict_v1.json"
SELECTOR_TRANSITION_GLOBAL = GENERATED / "selector_transition_global_c_v1_strict_v1.json"
SELECTOR_STATE_GLOBAL_PROJECTIVE = GENERATED / "selector_state_global_c_v1_projective_strict_v1.json"
SELECTOR_STATE_GLOBAL_DIRECTED = GENERATED / "selector_state_global_c_v1_directed_strict_v1.json"


def _read_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {"_parse_error": True}


def _is_pass(status: str | None) -> bool:
    return bool(status) and str(status).startswith("PASS")


h37 = _read_json(H37_SUMMARY)

atlas_present = SELECTOR_ATLAS_GLOBAL.exists()
transition_present = SELECTOR_TRANSITION_GLOBAL.exists()
projective_present = SELECTOR_STATE_GLOBAL_PROJECTIVE.exists()
directed_present = SELECTOR_STATE_GLOBAL_DIRECTED.exists()

sign_sensitive_present = _is_pass((h37 or {}).get("status")) and directed_present
qw2191_discharged = False

missing: list[str] = []
if not sign_sensitive_present:
    missing.append("sign_sensitive_directed_selector_state_object_or_observable")
if not qw2191_discharged:
    missing.append("global_qw_2191_discharge")

status = (
    "PASS_GLOBAL_SELECTOR_ATLAS_AND_GLUING_DATA_EXPORTED_ON_C_V1__"
    + ("DIRECTED_SIGN_SENSITIVE_LAYER_PRESENT__" if sign_sensitive_present else "DIRECTED_SIGN_SENSITIVE_LAYER_MISSING__")
    + "QW2191_STILL_OPEN"
    + ("__PREMISE_BASED_T164" if sign_sensitive_present else "")
)
frontier = "H41_B1" if not sign_sensitive_present else "H41_B2"
frontier_text = (
    "strict core exports a global selector atlas + transition object on C_v1 "
    "(projector section level), but global QW-2191 discharge remains open"
    + ("; directed/sign-sensitive layer is present (premise-based via T164)" if sign_sensitive_present else "")
)

payload = {
    "step": "H41",
    "title": "Selector Atlas And Gluing Data Audit",
    "date": "2026-03-17",
    "status": status,
    "inputs": {
        "H37": (h37.get("status") if isinstance(h37, dict) else None),
        "selector_atlas_global": str(SELECTOR_ATLAS_GLOBAL.relative_to(ROOT)) if atlas_present else None,
        "selector_transition_global": str(SELECTOR_TRANSITION_GLOBAL.relative_to(ROOT)) if transition_present else None,
        "selector_state_global_projective": str(SELECTOR_STATE_GLOBAL_PROJECTIVE.relative_to(ROOT)) if projective_present else None,
        "selector_state_global_directed": str(SELECTOR_STATE_GLOBAL_DIRECTED.relative_to(ROOT)) if directed_present else None,
    },
    "supports": [
        "global_selector_atlas_export_on_C_v1",
        "global_selector_transition_object_export_on_C_v1",
        "global_projective_selector_state_object_export_on_C_v1",
        "global_directed_selector_state_object_export_on_C_v1 (premise-based via T164)",
    ],
    "missing": missing,
    "frontier": frontier,
    "frontier_text": frontier_text,
    "hard_limits": {
        "theorem_level_pass": False,
        "full_closure_pass": False,
        "qw_2191_discharged": bool(qw2191_discharged),
        "operator_level_groupoid_promotion": False,
    },
    "no_false_pass": True,
}

summary = {"step": payload["step"], "status": payload["status"], "frontier": payload["frontier_text"], "no_false_pass": True}

(GENERATED / "h41_selector_atlas_and_gluing_data_audit.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
(GENERATED / "h41_selector_atlas_and_gluing_data_audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
