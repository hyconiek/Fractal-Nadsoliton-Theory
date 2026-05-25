#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2067_s1017_strict_same_scheme_certificate_transportability_targeting_audit.json"
MD = GEN / "p2067_s1017_strict_same_scheme_certificate_transportability_targeting_audit.md"

SCHEMA_VERSION = "p2067_s1017_v1"
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

    p2065 = load("p2065_s1015_strict_same_scheme_certificate_transportability_gap_audit.json")
    p2046 = load("p2046_s996_strict_same_scheme_adversarial_bucket_perturbation_audit.json")

    ready = p2065.get("result_kind") == "PASS_CERTIFICATE_TRANSPORTABILITY_GAP_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    gap_block = p2065.get("transportability_gap") or {}
    radius_star = float(gap_block.get("radius_star", 0.0))
    observed_amp = float(gap_block.get("observed_perturbation_amplitude_proxy", 0.0))
    transportability_gap = float(gap_block.get("transportability_gap", 0.0))

    eps_candidates: list[float] = []
    for key in ("epsilon_grid", "tested_epsilons", "adversarial_epsilons"):
        vals = p2046.get(key)
        if isinstance(vals, list):
            for value in vals:
                try:
                    eps_candidates.append(abs(float(value)))
                except Exception:
                    pass
    eps_max = max(eps_candidates, default=0.0)

    required_reduction_factor = (radius_star / observed_amp) if observed_amp > 0.0 else 1.0
    required_reduction_percent = max(0.0, (1.0 - required_reduction_factor) * 100.0)
    safe_next_amplitude_budget = max(0.0, radius_star)

    largest_eps_still_safe = min(eps_max, safe_next_amplitude_budget)
    currently_transportable = observed_amp <= radius_star

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2067",
        "stage_id": "S1017",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_CERTIFICATE_TRANSPORTABILITY_TARGETING_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready
            else "OPEN_CERTIFICATE_TRANSPORTABILITY_TARGETING_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2065_present": p2065.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2065_json_sha256": file_sha256(GEN / "p2065_s1015_strict_same_scheme_certificate_transportability_gap_audit.json"),
            "p2046_json_sha256": file_sha256(GEN / "p2046_s996_strict_same_scheme_adversarial_bucket_perturbation_audit.json"),
        },
        "transportability_targeting": {
            "radius_star": radius_star,
            "observed_perturbation_amplitude_proxy": observed_amp,
            "transportability_gap": transportability_gap,
            "currently_transportable": currently_transportable,
            "safe_next_amplitude_budget": safe_next_amplitude_budget,
            "required_reduction_factor": required_reduction_factor,
            "required_reduction_percent": required_reduction_percent,
            "largest_eps_candidate_still_safe": largest_eps_still_safe,
        },
        "c3_gate_update": {
            "C3_certificate_transportability_targeting_audit": "COMPUTED",
            "C3_transport_theorem": "OPEN",
            "C3_discharge_status": "NOT_DISCHARGED",
            "missing_for_c3_theorem": [
                "theorem-grade operator identity across full background family",
                "cross-background finite-part transport identity theorem",
                "global finite-part lock/cocycle theorem",
            ],
        },
        "gatekeeper_checks": {
            "preconditions_ready": ready,
            "gap_nonnegative": transportability_gap >= 0.0,
            "required_reduction_factor_in_0_1": 0.0 <= required_reduction_factor <= 1.0,
            "required_reduction_percent_nonnegative": required_reduction_percent >= 0.0,
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
                "# P2067 S1017: certificate transportability targeting audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- currently_transportable: `{currently_transportable}`",
                f"- required_reduction_factor: `{required_reduction_factor}`",
                f"- required_reduction_percent: `{required_reduction_percent}`",
                f"- safe_next_amplitude_budget: `{safe_next_amplitude_budget}`",
                "",
                "This stage turns the P2065 gap into a conservative targeting bound.",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
