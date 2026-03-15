#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

QW2122_JSON = REPO / "report_qw2122_psi_potential_diagonal_floor_gate.json"
ALPHA_GEO_JSON = GENERATED / "alpha_geo_strict_derived_v1.json"

OUT_JSON = (
    GENERATED
    / "ax25_strict_extension_lane_qw2122_shannon_directed_vacuum_and_self_couplings_value_provider_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "ax25_strict_extension_lane_qw2122_shannon_directed_vacuum_and_self_couplings_value_provider_packet_summary.json"
)


def is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def alpha_geo_numeric_from_export(obj: dict[str, Any]) -> float | None:
    value = obj.get("value")
    if is_number(value):
        return float(value)
    if isinstance(value, str):
        # Minimal parser for the currently exported form: "4 ln(2)".
        v = value.replace(" ", "")
        if v.lower() in {"4ln(2)", "ln(16)"}:
            return 4.0 * math.log(2.0)
    return None


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing_files: list[str] = []
    for p in (QW2122_JSON, ALPHA_GEO_JSON):
        if not p.exists():
            missing_files.append(str(p.relative_to(REPO)))

    if missing_files:
        summary = {
            "stage": "AX25",
            "status": "NOT_COMPUTABLE_MISSING_REQUIRED_FILES",
            "missing_files": missing_files,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    r2122 = load_json(QW2122_JSON)
    alpha_geo_export = load_json(ALPHA_GEO_JSON)

    lambda_psi_strict = (r2122.get("inputs") or {}).get("lambda_psi_strict")
    rho_star_sq = ((r2122.get("branch_results") or {}).get("broken_branch_strict") or {}).get("rho_star_sq")

    missing_fields: list[str] = []
    if not is_number(lambda_psi_strict):
        missing_fields.append("QW-2122.inputs.lambda_psi_strict (finite number)")
    if not is_number(rho_star_sq) or float(rho_star_sq) <= 0.0:
        missing_fields.append("QW-2122.branch_results.broken_branch_strict.rho_star_sq (positive finite number)")

    alpha_geo_value = alpha_geo_numeric_from_export(alpha_geo_export)
    if alpha_geo_value is None:
        missing_fields.append("alpha_geo_strict_derived_v1 numeric value (expected '4 ln(2)' export form)")

    if missing_fields:
        summary = {
            "stage": "AX25",
            "status": "NOT_COMPUTABLE_MISSING_REQUIRED_FIELDS",
            "missing_fields": missing_fields,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    lambda_psi_strict_f = float(lambda_psi_strict)
    rho_star_sq_f = float(rho_star_sq)
    assert alpha_geo_value is not None

    # Extension-only mapping premises (explicit):
    # 1) pick a marked-site + marked-direction coordinate on Z_12 via d_dir(0,i)=i (i=0..11),
    # 2) define weights w_i := exp(-alpha_geo * i) using strict-derived alpha_geo=4 ln 2,
    # 3) choose the vacuum squared amplitudes proportional to w_i:
    #       vpsi_i^2 := rho_*^2 * w_i / sum_j w_j
    #    so that Σ_i vpsi_i^2 = rho_*^2 exactly,
    # 4) keep the same uniform self-coupling mapping as AX24 (premise, not strict-derived):
    #       g4_i := 12*lambda_psi_strict,  g6_i := 0.
    weights = [math.exp(-alpha_geo_value * float(i)) for i in range(12)]
    sum_w = sum(weights)
    vpsi_sq = [rho_star_sq_f * (w / sum_w) for w in weights]
    vpsi = [math.sqrt(x) for x in vpsi_sq]

    g4_val = 12.0 * lambda_psi_strict_f
    g4 = [float(g4_val) for _ in range(12)]
    g6 = [0.0 for _ in range(12)]

    artifact = {
        "stage": "AX25",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": "STRICT_EXTENSION_ONLY_VALUE_PROVIDER_PACKET",
        "scope": "strict_extension_only",
        "goal": "provide_one_explicit_reproducible_candidate_value_provider_for_(vpsi,g4,g6)_using_strict_alpha_geo_directed_shannon_weights_to_enable_exploratory_P437/P434_execution_without_promoting_it_to_strict_core",
        "provenance": {
            "qw2122_source": str(QW2122_JSON.relative_to(REPO)),
            "qw2122_fields_used": {
                "lambda_psi_strict": "inputs.lambda_psi_strict",
                "rho_star_sq": "branch_results.broken_branch_strict.rho_star_sq",
            },
            "alpha_geo_source_upgrade": str(ALPHA_GEO_JSON.relative_to(REPO)),
            "alpha_geo_value_numeric": float(alpha_geo_value),
        },
        "explicit_mapping_premises": {
            "marked_site_and_direction": "use the directed coordinate i := d_dir(0,i)=i on Z_12 (marked site 0 + marked direction)",
            "directed_shannon_weights": "w_i := exp(-alpha_geo_strict_derived_v1 * i) with alpha_geo_strict_derived_v1 = 4 ln 2",
            "vacuum_from_weights": "vpsi_i^2 := rho_star_sq * w_i / sum_j w_j (so Σ_i vpsi_i^2 = rho_star_sq)",
            "uniform_self_couplings": "g4_i := 12*lambda_psi_strict, g6_i := 0 for i=0..11 (premise; not strict-derived)",
            "warning": "This packet is NOT a strict-derived discharge of T168; it adds explicit mapping premises to run evaluation harnesses.",
        },
        # Compatibility: P437 expects the numeric lists at top-level keys (vpsi,g4,g6).
        "vpsi": [float(x) for x in vpsi],
        "g4": g4,
        "g6": g6,
        "values": {
            "vpsi": [float(x) for x in vpsi],
            "g4": g4,
            "g6": g6,
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "AX25",
        "status": "PASS_EXTENSION_VALUE_PROVIDER_EXPORTED",
        "export": str(OUT_JSON.relative_to(REPO)),
        "notes": "Extension-only candidate values for vpsi/g4/g6 constructed from QW-2122 scalar vacuum closure + strict-derived alpha_geo via explicit marked-direction mapping premises.",
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

