#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2063_s1013_strict_same_scheme_cross_model_robustness_certificate_audit.json"
MD = GEN / "p2063_s1013_strict_same_scheme_cross_model_robustness_certificate_audit.md"

SCHEMA_VERSION = "p2063_s1013_v1"
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
    p2062 = load("p2062_s1012_strict_same_scheme_nonlinearity_stress_audit.json")

    ready = p2062.get("result_kind") == "PASS_NONLINEARITY_STRESS_AUDIT_WITH_BEST_SHARE_STABILITY_MAP__C3_STILL_OPEN"
    per_model = p2062.get("per_model_results") or []

    # gather model -> best share and available shares
    model_best: dict[str, float] = {}
    model_shares: dict[str, set[float]] = {}
    for mr in per_model:
        model = str(mr.get("model"))
        rows = mr.get("rows") or []
        shares = {float(r.get("share")) for r in rows}
        model_shares[model] = shares
        best = mr.get("best_row") or {}
        if "share" in best:
            model_best[model] = float(best["share"])

    models_all = sorted(model_best.keys())
    common_shares = sorted(set.intersection(*(model_shares[m] for m in models_all)) if models_all else set())

    dominant = (p2062.get("best_share_stability") or {}).get("dominant_best_share")
    dominant_share = float(dominant) if dominant is not None else None

    unanimous_models = [m for m, s in model_best.items() if dominant_share is not None and s == dominant_share]
    unanimous = len(unanimous_models) == len(models_all) and len(models_all) > 0

    # stability interval certificate: all common shares equal to dominant under unanimity -> point interval
    if unanimous and dominant_share is not None:
        stable_interval = {"low": dominant_share, "high": dominant_share}
        certificate_kind = "POINT_STABLE_SHARE_CERTIFICATE"
    else:
        stable_interval = {"low": None, "high": None}
        certificate_kind = "NO_UNANIMOUS_SHARE_CERTIFICATE"

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2063",
        "stage_id": "S1013",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_CROSS_MODEL_ROBUSTNESS_CERTIFICATE_AUDIT_WITH_TRACE__C3_STILL_OPEN"
            if ready and len(models_all) > 0
            else "OPEN_CROSS_MODEL_ROBUSTNESS_CERTIFICATE_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2062_present": p2062.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2062_json_sha256": file_sha256(GEN / "p2062_s1012_strict_same_scheme_nonlinearity_stress_audit.json"),
        },
        "cross_model_summary": {
            "models_considered": models_all,
            "model_best_share": model_best,
            "common_shares": common_shares,
            "dominant_share": dominant_share,
            "unanimous_models": unanimous_models,
            "unanimous": unanimous,
        },
        "robustness_certificate": {
            "certificate_kind": certificate_kind,
            "stable_share_interval": stable_interval,
            "stable_model_set": unanimous_models,
            "interpretation": "certificate is local/audit-level and not a theorem-grade C3 discharge",
        },
        "c3_gate_update": {
            "C3_cross_model_robustness_certificate_audit": "COMPUTED",
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
            "models_nonempty": len(models_all) > 0,
            "certificate_emitted": certificate_kind is not None,
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
                "# P2063 S1013: cross-model robustness certificate audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- Certificate kind: `{certificate_kind}`",
                f"- Dominant share: `{dominant_share}`",
                "",
                "Certificate remains audit-level; C3 not discharged.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
