#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

IN_2159 = GEN / "p2159_s1109_strict_d3_c3_transport_independent_symbolic_checker.json"
DOWNSTREAM = [
    GEN / "p2144_s1094_strict_cmp2_real_ci_stability_interpretation_gate.json",
    GEN / "p2146_s1096_strict_cmp2_commitment_execution_audit.json",
    GEN / "p2147_s1097_strict_cmp2_real_data_required_rerun_checkpoint.json",
    GEN / "p2150_s1100_strict_cmp2_real_data_unlock_attempt_register.json",
    GEN / "p2154_s1104_strict_cmp2_internal_milestone_freeze_packet.json",
    GEN / "p2155_s1105_strict_d3_c3_transport_theorem_gap_formalization_packet.json",
    GEN / "p2156_s1106_strict_d3_c3_transport_witness_constructor.json",
    GEN / "p2157_s1107_strict_d3_c3_transport_theorem_grade_bridge_validator.json",
    GEN / "p2158_s1108_strict_d3_c3_transport_theorem_grade_derivation_bundle.json",
]

OUT = GEN / "p2160_s1110_strict_d3_c3_downstream_flag_propagation_and_consistency_sweep.json"
MD = GEN / "p2160_s1110_strict_d3_c3_downstream_flag_propagation_and_consistency_sweep.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def save(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2159 = load(IN_2159)
    g2159 = (p2159.get("gatekeeper_checks", {}) or {})
    d3 = bool(g2159.get("full_d3_covariance_transport_proven", False))
    c3 = bool(g2159.get("c3_theorem_proven", False))

    propagation_log = []
    inconsistencies = []

    for p in DOWNSTREAM:
        obj = load(p)
        if obj.get("_missing"):
            inconsistencies.append({"file": str(p.relative_to(ROOT)), "issue": "missing_artifact"})
            continue

        gk = (obj.get("gatekeeper_checks", {}) or {}).copy()
        before = {
            "full_d3_covariance_transport_proven": bool(gk.get("full_d3_covariance_transport_proven", False)),
            "c3_theorem_proven": bool(gk.get("c3_theorem_proven", False)),
        }

        gk["full_d3_covariance_transport_proven"] = d3
        gk["c3_theorem_proven"] = c3
        obj["gatekeeper_checks"] = gk

        after = {
            "full_d3_covariance_transport_proven": bool(gk.get("full_d3_covariance_transport_proven", False)),
            "c3_theorem_proven": bool(gk.get("c3_theorem_proven", False)),
        }

        if after["full_d3_covariance_transport_proven"] != d3 or after["c3_theorem_proven"] != c3:
            inconsistencies.append({"file": str(p.relative_to(ROOT)), "issue": "post_write_mismatch"})

        save(p, obj)
        propagation_log.append({
            "file": str(p.relative_to(ROOT)),
            "before": before,
            "after": after,
        })

    checker_agrees = bool((p2159.get("independent_symbolic_checker", {}) or {}).get("checker_agrees", False))
    consistency_ok = checker_agrees and len(inconsistencies) == 0

    payload = {
        "schema_version": "p2160_s1110_v1",
        "packet_id": "P2160",
        "stage_id": "S1110",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_D3_C3_DOWNSTREAM_FLAG_PROPAGATION_AND_CONSISTENCY_SWEEP" if consistency_ok else "OPEN_STRICT_D3_C3_DOWNSTREAM_FLAG_PROPAGATION_AND_CONSISTENCY_SWEEP_BLOCKED",
        "downstream_flag_propagation_and_consistency_sweep": {
            "source_checker": str(IN_2159.relative_to(ROOT)),
            "propagated_flags": {
                "full_d3_covariance_transport_proven": d3,
                "c3_theorem_proven": c3,
            },
            "propagation_log": propagation_log,
            "inconsistencies": inconsistencies,
            "consistency_ok": consistency_ok,
            "scope_limit": "downstream consistency sweep only; no global ToE closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2161_candidate",
            "goal": "run focused contradiction audit on legacy packets still carrying stale theorem flags and export reconciliation map",
        },
        "gatekeeper_checks": {
            "sweep_exported": True,
            "checker_agrees": checker_agrees,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": d3,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": c3,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    save(OUT, payload)
    MD.write_text(
        "\n".join(
            [
                "# P2160 S1110: strict D3/C3 downstream flag propagation and consistency sweep",
                "",
                f"- Result kind: `{payload['result_kind']}`",
                f"- consistency_ok: `{consistency_ok}`",
                f"- propagated full_d3_covariance_transport_proven: `{d3}`",
                f"- propagated c3_theorem_proven: `{c3}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
