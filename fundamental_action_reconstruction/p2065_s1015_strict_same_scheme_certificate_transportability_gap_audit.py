#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2065_s1015_strict_same_scheme_certificate_transportability_gap_audit.json"
MD = GEN / "p2065_s1015_strict_same_scheme_certificate_transportability_gap_audit.md"

SCHEMA_VERSION = "p2065_s1015_v1"
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
    p2064 = load("p2064_s1014_strict_same_scheme_certificate_perturbation_radius_audit.json")
    p2046 = load("p2046_s996_strict_same_scheme_adversarial_bucket_perturbation_audit.json")
    p2058 = load("p2058_s1008_strict_same_scheme_policy_regime_switch_audit.json")

    ready = p2064.get("result_kind") == "PASS_CERTIFICATE_PERTURBATION_RADIUS_AUDIT_WITH_TRACE__C3_STILL_OPEN"

    radius_star = float(((p2064.get("radius_scan") or {}).get("radius_star", 0.0)))

    # Real perturbation amplitude proxies from earlier audits.
    # 1) adversarial epsilon candidates from P2046 if present
    eps_candidates = []
    if isinstance(p2046, dict):
        for k in ("epsilon_grid", "tested_epsilons", "adversarial_epsilons"):
            vals = p2046.get(k)
            if isinstance(vals, list):
                for v in vals:
                    try:
                        eps_candidates.append(abs(float(v)))
                    except Exception:
                        pass

    # 2) regime parameter scale from P2058 grid (tau spread)
    tau_grid = ((p2058.get("regime_grid") or {}).get("softmax_tau") or []) if isinstance(p2058, dict) else []
    tau_amp = max([abs(float(t)) for t in tau_grid], default=0.0)

    observed_amp = max(eps_candidates + [tau_amp]) if (eps_candidates or tau_amp > 0) else 0.0

    transport_gap = max(0.0, observed_amp - radius_star)
    coverage_ratio = (radius_star / observed_amp) if observed_amp > 0 else 1.0

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2065",
        "stage_id": "S1015",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_CERTIFICATE_TRANSPORTABILITY_GAP_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready
            else "OPEN_CERTIFICATE_TRANSPORTABILITY_GAP_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2064_present": p2064.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2064_json_sha256": file_sha256(GEN / "p2064_s1014_strict_same_scheme_certificate_perturbation_radius_audit.json"),
            "p2046_json_sha256": file_sha256(GEN / "p2046_s996_strict_same_scheme_adversarial_bucket_perturbation_audit.json"),
            "p2058_json_sha256": file_sha256(GEN / "p2058_s1008_strict_same_scheme_policy_regime_switch_audit.json"),
        },
        "transportability_gap": {
            "radius_star": radius_star,
            "observed_perturbation_amplitude_proxy": observed_amp,
            "transportability_gap": transport_gap,
            "coverage_ratio": coverage_ratio,
        },
        "source_proxies": {
            "p2046_epsilon_candidates_count": len(eps_candidates),
            "p2058_tau_grid": tau_grid,
            "p2058_tau_amplitude_proxy": tau_amp,
        },
        "c3_gate_update": {
            "C3_certificate_transportability_gap_audit": "COMPUTED",
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
            "gap_nonnegative": transport_gap >= 0.0,
            "coverage_ratio_nonnegative": coverage_ratio >= 0.0,
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
                "# P2065 S1015: certificate transportability gap audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- radius_star: `{radius_star}`",
                f"- observed amplitude proxy: `{observed_amp}`",
                f"- transportability gap: `{transport_gap}`",
                "",
                "Gap quantifies what is still missing for wider certificate usability.",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
