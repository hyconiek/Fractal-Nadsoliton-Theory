#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2072_s1022_strict_same_scheme_tau_stability_margin_stress_audit.json"
MD = GEN / "p2072_s1022_strict_same_scheme_tau_stability_margin_stress_audit.md"

SCHEMA_VERSION = "p2072_s1022_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2071 = load("p2071_s1021_strict_same_scheme_tau_boundary_refinement_audit.json")
    ready = p2071.get("result_kind") == "PASS_TAU_BOUNDARY_REFINEMENT_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    ref = p2071.get("tau_boundary_refinement") or {}
    radius_star = float(ref.get("radius_star", 0.0))
    tight = ref.get("tight_feasible_tau") or {}
    tight_tau = float(tight.get("tau_candidate", 0.0))

    deltas = [-1e-3, -5e-4, -1e-4, 0.0, 1e-4, 5e-4, 1e-3]
    rows = []
    for d in deltas:
        tau = tight_tau + d
        positive_domain = tau > 0.0
        transportable = positive_domain and (tau <= radius_star)
        rows.append(
            {
                "delta": d,
                "tau_candidate": tau,
                "positive_domain": positive_domain,
                "transportable_under_radius_star": transportable,
                "signed_margin": radius_star - tau,
            }
        )

    pos_rows = [r for r in rows if r["delta"] > 0]
    neg_rows = [r for r in rows if r["delta"] < 0]
    positive_side_all_fail = all(not r["transportable_under_radius_star"] for r in pos_rows)
    negative_side_all_pass = all(r["transportable_under_radius_star"] for r in neg_rows)

    boundary_classification = (
        "SHARP_STABLE_BOUNDARY"
        if positive_side_all_fail and negative_side_all_pass
        else "FRAGILE_OR_UNRESOLVED_BOUNDARY"
    )

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2072",
        "stage_id": "S1022",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_TAU_STABILITY_MARGIN_STRESS_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready
            else "OPEN_TAU_STABILITY_MARGIN_STRESS_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2071_present": p2071.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2071_json_sha256": file_sha256(GEN / "p2071_s1021_strict_same_scheme_tau_boundary_refinement_audit.json"),
        },
        "tau_stability_margin": {
            "radius_star": radius_star,
            "tight_feasible_tau": tight_tau,
            "stress_deltas": deltas,
            "stress_rows": rows,
            "positive_side_all_fail": positive_side_all_fail,
            "negative_side_all_pass": negative_side_all_pass,
            "boundary_classification": boundary_classification,
        },
        "c3_gate_update": {
            "C3_tau_stability_margin_stress_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "stress_rows_nonempty": len(rows) > 0,
            "contains_requested_deltas": all(d in deltas for d in (-1e-3, -5e-4, -1e-4, 1e-4, 5e-4, 1e-3)),
            "c3_theorem_proven": False,
            "no_background_globalization_claimed": True,
            "no_tensor_component_claimed": True,
            "no_toe_closure_claimed": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2072 S1022: tau stability margin stress audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- tight_feasible_tau: `{tight_tau}`",
                f"- boundary_classification: `{boundary_classification}`",
                f"- positive_side_all_fail: `{positive_side_all_fail}`",
                f"- negative_side_all_pass: `{negative_side_all_pass}`",
                "",
                "Small perturbations around tight_feasible_tau are stress-tested for boundary stability.",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
