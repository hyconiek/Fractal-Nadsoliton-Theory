#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

QW2122_JSON = REPO / "report_qw2122_psi_potential_diagonal_floor_gate.json"

R_ORD_JSON = GENERATED / "r_ord_z12_v1_reference_distribution.json"

IN_P437 = GENERATED / "p447_input_vpsi_g4_g6_element_order_reference_candidate.json"
OUT_SUMMARY = GENERATED / "p447_current_strict_t169_element_order_reference_lift_candidate_pipeline_audit_probe_summary.json"

OUT_P437 = GENERATED / "p447_p437_out.json"
OUT_P437_SUMMARY = GENERATED / "p447_p437_out_summary.json"

IN_P434 = GENERATED / "p447_p434_input_sigma_values_from_p437.json"
OUT_P434 = GENERATED / "p447_p434_out.json"
OUT_P434_SUMMARY = GENERATED / "p447_p434_out_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def zn_element_order(n: int, k: int) -> int:
    kk = k % n
    if kk == 0:
        return 1
    return n // math.gcd(kk, n)


def normalize_positive_weights(w: list[float]) -> list[float]:
    if any((not is_number(x)) or float(x) <= 0.0 for x in w):
        raise ValueError("weights must be finite and strictly positive")
    z = float(sum(float(x) for x in w))
    if z <= 0.0:
        raise ValueError("weights must sum to positive")
    return [float(x) / z for x in w]


def run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, check=False)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing_files: list[str] = []
    for p in [QW2122_JSON, R_ORD_JSON]:
        if not p.exists():
            missing_files.append(str(p.relative_to(REPO) if p.is_absolute() else p))

    if missing_files:
        summary = {
            "stage": "P447",
            "status": "NOT_COMPUTABLE_MISSING_REQUIRED_FILES",
            "missing_files": missing_files,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    r2122 = load_json(QW2122_JSON)
    rho_star_sq = ((r2122.get("branch_results") or {}).get("broken_branch_strict") or {}).get("rho_star_sq")
    lambda_psi_strict = (r2122.get("inputs") or {}).get("lambda_psi_strict")

    missing_fields: list[str] = []
    if not is_number(rho_star_sq):
        missing_fields.append("QW-2122.branch_results.broken_branch_strict.rho_star_sq (finite number)")
    if not is_number(lambda_psi_strict):
        missing_fields.append("QW-2122.inputs.lambda_psi_strict (finite number)")

    if missing_fields:
        summary = {
            "stage": "P447",
            "status": "NOT_COMPUTABLE_MISSING_REQUIRED_FIELDS",
            "missing_fields": missing_fields,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    rho_star_sq_f = float(rho_star_sq)
    lambda_psi_strict_f = float(lambda_psi_strict)

    # Use the strict-derived alpha_geo value (symbolic "4 ln 2") in numeric form.
    alpha_geo = 4.0 * math.log(2.0)

    # Element-order reference on Z_12 (direction-free by N479 / exported by F446).
    n = 12
    orders = [zn_element_order(n, i) for i in range(n)]
    weights = [math.exp(-alpha_geo * float(o)) for o in orders]
    rprob = normalize_positive_weights(weights)

    # Explicit non-strict mapping premises (diagnostic only):
    # vpsi_i^2 := rho_star_sq * rprob_i,  g4_i := 12*lambda_psi_strict,  g6_i := 0.
    vpsi = [math.sqrt(rho_star_sq_f * float(r)) for r in rprob]
    g4_val = 12.0 * lambda_psi_strict_f
    g4 = [float(g4_val) for _ in range(n)]
    g6 = [0.0 for _ in range(n)]

    IN_P437.write_text(
        json.dumps(
            {
                "stage": "P447",
                "generated_utc": datetime.now(timezone.utc).isoformat(),
                "status": "CANDIDATE_INPUT_FOR_P437_DIAGNOSTIC_ONLY",
                "scope": "probe_only_non_strict_mapping_premises",
                "premises": {
                    "vpsi_magnitudes": "vpsi_i^2 := rho_star_sq * r_ord(i) (element-order reference; no marked direction)",
                    "uniform_self_couplings": "g4_i := 12*lambda_psi_strict, g6_i := 0",
                },
                "inputs": {
                    "QW-2122": str(QW2122_JSON.relative_to(REPO)),
                    "F446_reference": str(R_ORD_JSON.relative_to(REPO)),
                    "alpha_geo_used_numeric": alpha_geo,
                },
                "derived_reference": {"ord_Z12": orders, "reference_prob": rprob},
                "vpsi": [float(x) for x in vpsi],
                "g4": g4,
                "g6": g6,
                "no_false_pass": True,
            },
            indent=2,
            ensure_ascii=True,
        )
        + "\n",
        encoding="ascii",
    )

    # Run P437 (N477 evaluation harness).
    proc_p437 = run(
        [
            "python3",
            "fundamental_action_reconstruction/p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe.py",
            "--in",
            str(IN_P437.relative_to(REPO)),
            "--out-json",
            str(OUT_P437.relative_to(REPO)),
            "--out-summary",
            str(OUT_P437_SUMMARY.relative_to(REPO)),
        ]
    )

    if proc_p437.returncode != 0 or not OUT_P437.exists():
        summary = {
            "stage": "P447",
            "status": "FAILED_TO_RUN_P437",
            "p437_returncode": proc_p437.returncode,
            "p437_stdout": proc_p437.stdout[-2000:],
            "p437_stderr": proc_p437.stderr[-2000:],
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p437_out = load_json(OUT_P437)
    sigmas = (((p437_out.get("computed") or {}).get("sigma_opposite_pair_sums")) or {}) if isinstance(p437_out, dict) else {}

    required = [
        "Sigma_psi0_psi6",
        "Sigma_psi1_psi7",
        "Sigma_psi2_psi8",
        "Sigma_psi3_psi9",
        "Sigma_psi4_psi10",
        "Sigma_psi5_psi11",
    ]
    missing_sigmas = [k for k in required if not is_number((sigmas or {}).get(k))]
    if missing_sigmas:
        summary = {
            "stage": "P447",
            "status": "P437_OUTPUT_MISSING_SIGMAS",
            "missing_sigma_keys": missing_sigmas,
            "p437_out": str(OUT_P437.relative_to(REPO)),
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    IN_P434.write_text(
        json.dumps({k: float(sigmas[k]) for k in required}, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )

    # Run P434 evaluation.
    proc_p434 = run(
        [
            "python3",
            "fundamental_action_reconstruction/p434_current_strict_canonical_local_diagonal_mode2_defect_value_instantiation_evaluation_probe.py",
            "--in",
            str(IN_P434.relative_to(REPO)),
            "--out-json",
            str(OUT_P434.relative_to(REPO)),
            "--out-summary",
            str(OUT_P434_SUMMARY.relative_to(REPO)),
        ]
    )

    if proc_p434.returncode != 0 or not OUT_P434.exists():
        summary = {
            "stage": "P447",
            "status": "FAILED_TO_RUN_P434",
            "p434_returncode": proc_p434.returncode,
            "p434_stdout": proc_p434.stdout[-2000:],
            "p434_stderr": proc_p434.stderr[-2000:],
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p434_out = load_json(OUT_P434)
    result = (p434_out.get("result") or {}) if isinstance(p434_out, dict) else {}

    summary = {
        "stage": "P447",
        "status": "EXECUTED_DIAGNOSTIC_PIPELINE_NO_FALSE_PASS",
        "inputs": {
            "QW-2122": str(QW2122_JSON.relative_to(REPO)),
            "F446_reference": str(R_ORD_JSON.relative_to(REPO)),
            "P437_in": str(IN_P437.relative_to(REPO)),
        },
        "outputs": {
            "P437_out": str(OUT_P437.relative_to(REPO)),
            "P434_out": str(OUT_P434.relative_to(REPO)),
        },
        "p434_result": {
            "cuts_O2_on_pair1_by_N466": result.get("cuts_O2_on_pair1_by_N466"),
            "F2_abs": (result.get("F2") or {}).get("abs") if isinstance(result.get("F2"), dict) else None,
            "theta_star_by_N468": result.get("theta_star_by_N468"),
        },
        "notes": "This is probe-level only: it uses explicit non-strict mapping premises to drive P437/P434.",
        "no_false_pass": True,
    }

    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

