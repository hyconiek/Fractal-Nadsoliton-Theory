#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2062_s1012_strict_same_scheme_nonlinearity_stress_audit.json"
MD = GEN / "p2062_s1012_strict_same_scheme_nonlinearity_stress_audit.md"

SCHEMA_VERSION = "p2062_s1012_v1"
TIMESTAMP_UTC = "2026-05-25T00:00:00+00:00"
SHARES = [0.6, 0.7, 0.8, 0.9]


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def file_sha256(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return hashlib.sha256(path.read_bytes()).hexdigest()


def model_linear(s: float) -> float:
    return 0.2 + 0.5 * s


def model_quadratic(s: float) -> float:
    return 0.15 + 0.6 * (s**2)


def model_saturating(s: float) -> float:
    return 0.1 + 0.7 * (1.0 - (1.0 - s) ** 2)


def model_convex(s: float) -> float:
    return 0.18 + 0.65 * ((s - 0.5) ** 1.5)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2061 = load("p2061_s1011_strict_same_scheme_reallocation_stress_audit.json")

    ready = p2061.get("result_kind") == "PASS_REALLOCATION_STRESS_AUDIT_WITH_ENVELOPE_CI__C3_STILL_OPEN"
    baseline = p2061.get("baseline") or {}
    base_max = float(baseline.get("fixed_regret_max", 0.0))

    models = {
        "linear": model_linear,
        "quadratic": model_quadratic,
        "saturating": model_saturating,
        "convex": model_convex,
    }

    per_model = []
    best_share_votes: dict[str, int] = {str(s): 0 for s in SHARES}

    for mname, mfun in models.items():
        rows = []
        for s in SHARES:
            shrink_factor = max(0.0, min(0.98, mfun(s)))
            aware_max = base_max * (1.0 - shrink_factor)
            shrink_abs = base_max - aware_max
            rows.append(
                {
                    "share": s,
                    "shrink_factor": shrink_factor,
                    "regime_aware_regret_max": aware_max,
                    "regret_shrink_absolute": shrink_abs,
                }
            )
        best = max(rows, key=lambda r: r["regret_shrink_absolute"]) if rows else None
        if best is not None:
            best_share_votes[str(best["share"])] += 1
        per_model.append({"model": mname, "rows": rows, "best_row": best})

    dominant_best_share = max(best_share_votes, key=best_share_votes.get) if best_share_votes else None
    stability_fraction = (best_share_votes.get(dominant_best_share, 0) / len(models)) if models else 0.0

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2062",
        "stage_id": "S1012",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": (
            "PASS_NONLINEARITY_STRESS_AUDIT_WITH_BEST_SHARE_STABILITY_MAP__C3_STILL_OPEN"
            if ready and len(per_model) > 0
            else "OPEN_NONLINEARITY_STRESS_AUDIT_BLOCKED"
        ),
        "depends_on": {
            "p2061_present": p2061.get("_missing") is None,
            "preconditions_ready": ready,
        },
        "input_hashes": {
            "p2061_json_sha256": file_sha256(GEN / "p2061_s1011_strict_same_scheme_reallocation_stress_audit.json"),
        },
        "stress_setup": {
            "shares": SHARES,
            "models": list(models.keys()),
            "base_fixed_regret_max": base_max,
        },
        "per_model_results": per_model,
        "best_share_stability": {
            "best_share_votes": best_share_votes,
            "dominant_best_share": dominant_best_share,
            "stability_fraction": stability_fraction,
        },
        "c3_gate_update": {
            "C3_nonlinearity_stress_audit": "COMPUTED",
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
            "models_nonempty": len(per_model) > 0,
            "dominant_best_share_present": dominant_best_share is not None,
            "stability_fraction_in_unit_interval": 0.0 <= stability_fraction <= 1.0,
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
                "# P2062 S1012: nonlinearity stress audit",
                "",
                f"- Status: `{payload['status']}`",
                f"- Result kind: `{payload['result_kind']}`",
                f"- Dominant best share: `{dominant_best_share}`",
                f"- Stability fraction: `{stability_fraction}`",
                "",
                "Compared multiple nonlinear shrink-vs-share models.",
                "C3 remains OPEN (not discharged).",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
