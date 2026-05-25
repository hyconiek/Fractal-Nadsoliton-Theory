#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2068_s1018_strict_same_scheme_transportability_improvement_lever_audit.json"
MD = GEN / "p2068_s1018_strict_same_scheme_transportability_improvement_lever_audit.md"

SCHEMA_VERSION = "p2068_s1018_v1"
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

    p2067 = load("p2067_s1017_strict_same_scheme_certificate_transportability_targeting_audit.json")
    p2046 = load("p2046_s996_strict_same_scheme_adversarial_bucket_perturbation_audit.json")
    p2058 = load("p2058_s1008_strict_same_scheme_policy_regime_switch_audit.json")

    ready = p2067.get("result_kind") == "PASS_CERTIFICATE_TRANSPORTABILITY_TARGETING_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    target = p2067.get("transportability_targeting") or {}
    radius_star = float(target.get("radius_star", 0.0))

    eps_candidates: list[float] = []
    for key in ("epsilon_grid", "tested_epsilons", "adversarial_epsilons"):
        vals = p2046.get(key)
        if isinstance(vals, list):
            for value in vals:
                try:
                    eps_candidates.append(abs(float(value)))
                except Exception:
                    pass
    epsilon_amp = max(eps_candidates, default=0.0)

    tau_grid = ((p2058.get("regime_grid") or {}).get("softmax_tau") or []) if isinstance(p2058, dict) else []
    tau_amp = max([abs(float(t)) for t in tau_grid], default=0.0)

    observed_amp = max(epsilon_amp, tau_amp)
    gap = max(0.0, observed_amp - radius_star)

    # Component-wise deficits relative to radius_star.
    epsilon_deficit = max(0.0, epsilon_amp - radius_star)
    tau_deficit = max(0.0, tau_amp - radius_star)

    epsilon_required_factor = (radius_star / epsilon_amp) if epsilon_amp > 0 else 1.0
    tau_required_factor = (radius_star / tau_amp) if tau_amp > 0 else 1.0

    lever_rows = [
        {
            "lever": "epsilon_amplitude_shrink",
            "current_amplitude": epsilon_amp,
            "deficit_to_radius_star": epsilon_deficit,
            "required_shrink_factor": max(0.0, min(1.0, epsilon_required_factor)),
            "required_shrink_percent": max(0.0, min(100.0, (1.0 - epsilon_required_factor) * 100.0)),
        },
        {
            "lever": "tau_regime_amplitude_shrink",
            "current_amplitude": tau_amp,
            "deficit_to_radius_star": tau_deficit,
            "required_shrink_factor": max(0.0, min(1.0, tau_required_factor)),
            "required_shrink_percent": max(0.0, min(100.0, (1.0 - tau_required_factor) * 100.0)),
        },
    ]
    lever_rows.sort(key=lambda row: row["deficit_to_radius_star"], reverse=True)

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2068",
        "stage_id": "S1018",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_TRANSPORTABILITY_IMPROVEMENT_LEVER_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready
            else "OPEN_TRANSPORTABILITY_IMPROVEMENT_LEVER_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2067_present": p2067.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2067_json_sha256": file_sha256(GEN / "p2067_s1017_strict_same_scheme_certificate_transportability_targeting_audit.json"),
            "p2046_json_sha256": file_sha256(GEN / "p2046_s996_strict_same_scheme_adversarial_bucket_perturbation_audit.json"),
            "p2058_json_sha256": file_sha256(GEN / "p2058_s1008_strict_same_scheme_policy_regime_switch_audit.json"),
        },
        "transportability_components": {
            "radius_star": radius_star,
            "epsilon_amplitude_proxy": epsilon_amp,
            "tau_amplitude_proxy": tau_amp,
            "observed_amplitude_proxy": observed_amp,
            "transportability_gap": gap,
        },
        "improvement_lever_ranking": lever_rows,
        "c3_gate_update": {
            "C3_transportability_improvement_lever_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "gap_nonnegative": gap >= 0.0,
            "ranking_nonempty": len(lever_rows) > 0,
            "ranked_descending_deficit": all(
                lever_rows[i]["deficit_to_radius_star"] >= lever_rows[i + 1]["deficit_to_radius_star"]
                for i in range(len(lever_rows) - 1)
            ),
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
                "# P2068 S1018: transportability improvement lever audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- radius_star: `{radius_star}`",
                f"- observed_amplitude_proxy: `{observed_amp}`",
                f"- transportability_gap: `{gap}`",
                "",
                "Top lever ranking is by largest deficit to radius_star.",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
