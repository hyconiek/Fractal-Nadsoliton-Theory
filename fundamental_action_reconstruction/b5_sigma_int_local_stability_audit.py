#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

AS_OF = "2026-03-15"

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
GAUGE_QUOTIENT_SAFETY_WITNESS = GENERATED / "sigma_int_gauge_quotient_safety_witness_v1.json"


def load_if_exists(path: Path) -> dict | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


gauge_witness = load_if_exists(GAUGE_QUOTIENT_SAFETY_WITNESS)
gauge_quotient_safety_exported = bool(
    gauge_witness is not None and gauge_witness.get("object") == "sigma_int_gauge_quotient_safety_witness_v1"
)

status = "B5_EXECUTED_LOCAL_DEFORMATION_STABILITY_SUPPORTED_GAUGE_QUOTIENT_DISCHARGE_PENDING_NO_FALSE_PASS"
full_gauge_quotient_status = "open"
if gauge_quotient_safety_exported:
    status = "B5_UPDATED_LOCAL_DEFORMATION_STABILITY_SUPPORTED_GAUGE_QUOTIENT_SAFETY_WITNESS_EXPORTED_NO_FALSE_PASS"
    full_gauge_quotient_status = "discharged_on_declared_domain"

summary = {
    "program_id": "fundamental_action_reconstruction",
    "as_of": AS_OF,
    "status": status,
    "source_policy": {
        "strict_admissible_support": [
            "QW-2206",
            "F308/N419",
            "A2",
            "A3",
            "B4",
            "A10",
        ],
        "hybrid_support": [
            "QW-1622",
        ],
    },
    "anti_overclaim": {
        "b3_o2_discharged_claim": False,
        "full_gauge_safety_claim": gauge_quotient_safety_exported,
        "axiom_free_uniqueness_closed_claim": False,
    },
    "b5": {
        "goal": "Assess whether sigma_int_candidate survives local topological deformations and whether gauge-choice dependence is already ruled out.",
        "findings": [
            {
                "object": "local deformation stability in fixed topological sector",
                "status": "supported_partial",
                "basis": ["QW-2206", "QW-1622"],
            },
            {
                "object": "independence from simple parametrization choice",
                "status": "supported_partial",
                "basis": ["B4"],
            },
            {
                "object": "full gauge quotient safety",
                "status": full_gauge_quotient_status,
                "basis": (["F308/N419"] if gauge_quotient_safety_exported else ["A2", "A3"]),
                "witness_artifact": (
                    "generated/sigma_int_gauge_quotient_safety_witness_v1.json" if gauge_quotient_safety_exported else None
                ),
            },
            {
                "object": "B3_O2 discharge",
                "status": (
                    "discharged_on_declared_domain"
                    if gauge_quotient_safety_exported
                    else "partial_local_support_only"
                ),
            },
        ],
        "obligation_status": {
            "B3_O1": "candidate_identified",
            "B3_O2": (
                "discharged_on_declared_domain"
                if gauge_quotient_safety_exported
                else "partial_local_support_only"
            ),
            "B3_O3": "open",
            "B3_O4": "open",
            "B3_O5": "pending_after_O3_to_O4",
        },
        "next_step": "B6",
    },
}

out = ROOT / "generated" / "b5_sigma_int_local_stability_audit_summary.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
print(f"wrote {out}")
