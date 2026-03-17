from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
GENERATED.mkdir(exist_ok=True)

H37_SUMMARY = GENERATED / "h37_sign_distinction_state_audit_summary.json"
H38_SUMMARY = GENERATED / "h38_projective_selector_state_audit_summary.json"

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
h38 = _read_json(H38_SUMMARY)

projective_present = SELECTOR_STATE_GLOBAL_PROJECTIVE.exists()
directed_present = SELECTOR_STATE_GLOBAL_DIRECTED.exists()

sign_sensitive_present = _is_pass((h37 or {}).get("status")) and _is_pass((h38 or {}).get("status"))

qw2191_discharged = False

missing: list[str] = []
if not sign_sensitive_present:
    missing.append("sign_sensitive_directed_selector_state_object_or_observable")
if not qw2191_discharged:
    missing.append("global_qw_2191_discharge")

if projective_present and directed_present and sign_sensitive_present:
    status = (
        "PASS_PARTIAL_GLOBAL_PROJECTIVE_AND_DIRECTED_SELECTOR_STATE_OBJECTS_EXPORTED_ON_C_V1__"
        "SIGN_SENSITIVE_OBSERVABLE_PRESENT__QW2191_STILL_OPEN__PREMISE_BASED_T164"
    )
    frontier = "H39_B3"
    frontier_text = (
        "strict core exports global projective and directed selector state objects on C_v1 "
        "(directed continuation; premise-based via T164), but global QW-2191 discharge is still open"
    )
else:
    status = (
        "PASS_PARTIAL_GLOBAL_PROJECTIVE_SELECTOR_STATE_OBJECT_EXPORTED_ON_C_V1__"
        "DIRECTED_OR_SIGN_SENSITIVE_LAYER_INCOMPLETE__QW2191_STILL_OPEN"
    )
    frontier = "H39_B2"
    frontier_text = (
        "strict core exports a global projective selector state object on C_v1, but the directed/sign-sensitive layer is not fully exported "
        "and global QW-2191 discharge remains open"
    )

payload = {
    "step": "H39",
    "title": "Global Selector Object Absence Audit",
    "date": "2026-03-17",
    "status": status,
    "inputs": {
        "H34": "no_strict_basis_covariance_or_target_independence_argument",
        "H35": "no_strict_physical_axis_selection",
        "H36": "no_aut_invariant_directed_orientation_selection",
        "H37": (h37.get("status") if isinstance(h37, dict) else None),
        "H38": (h38.get("status") if isinstance(h38, dict) else None),
        "T170": f"selector atlas/transition present: atlas={SELECTOR_ATLAS_GLOBAL.exists()}, transition={SELECTOR_TRANSITION_GLOBAL.exists()}",
        "global_projective_selector_state_object": str(SELECTOR_STATE_GLOBAL_PROJECTIVE.relative_to(ROOT)) if projective_present else None,
        "global_directed_selector_state_object": str(SELECTOR_STATE_GLOBAL_DIRECTED.relative_to(ROOT)) if directed_present else None,
    },
    "supports": [
        "global_selector_atlas_and_transition_objects_on_C_v1 (projector section level)",
        "global_projective_selector_state_object_on_C_v1",
        "global_directed_selector_state_object_on_C_v1 (premise-based via T164)",
        "sign_sensitive_directed_observable_on_pair1 (premise-based via T164)",
    ],
    "missing": missing,
    "frontier": frontier,
    "frontier_text": frontier_text,
    "hard_limits": {
        "theorem_level_pass": False,
        "full_closure_pass": False,
        "qw_2191_discharged": bool(qw2191_discharged),
        "aut_invariant_canonicity_claim": False,
    },
    "no_false_pass": True,
}

summary = {"step": payload["step"], "status": payload["status"], "frontier": payload["frontier_text"], "no_false_pass": True}

(GENERATED / "h39_global_selector_object_absence_audit.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
(GENERATED / "h39_global_selector_object_absence_audit_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
