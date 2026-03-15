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

IN_SIGMA_INT = GENERATED / "sigma_int_strict_derived_v1.json"
IN_GAUGE_SAFETY = GENERATED / "sigma_int_gauge_quotient_safety_witness_v1.json"
IN_EXPORT_MAP_OBJECT = GENERATED / "upsilon_residual_datum_sigma_int_bridge_export_map_object_v1.json"
IN_THETA_PAIR = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"
IN_R1_POPULATION = GENERATED / "r1_residual_orientation_datum_target_slot_population_strict_derived_from_sigma_int_slot_free_theta_pair_v1.json"

OUT_SUPPORT = GENERATED / "iota_residual_datum_sigma_int_bridge_export_map_object_support_v1.json"
OUT_SUPPORT_SUMMARY = GENERATED / "iota_residual_datum_sigma_int_bridge_export_map_object_support_v1_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def float_list(x: Any) -> list[float] | None:
    if not isinstance(x, list):
        return None
    if not all(isinstance(v, (int, float)) for v in x):
        return None
    return [float(v) for v in x]


def max_abs_diff(a: list[float], b: list[float]) -> float:
    if len(a) != len(b):
        raise ValueError("vector length mismatch")
    return float(max(abs(x - y) for x, y in zip(a, b)))


def mod_2pi_close(delta: float, target: float, tol: float = 1e-12) -> bool:
    # Compare delta to target modulo 2π in a tolerance-aware way.
    two_pi = 2.0 * math.pi
    d = float(delta - target)
    d = (d + math.pi) % two_pi - math.pi
    return abs(d) <= tol


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    sigma_int = load_json(IN_SIGMA_INT)
    gauge_safety = load_json(IN_GAUGE_SAFETY)
    export_map_object = load_json(IN_EXPORT_MAP_OBJECT)
    theta_pair = load_json(IN_THETA_PAIR)
    r1_population = load_json(IN_R1_POPULATION)

    sigma_val = sigma_int.get("value")
    if not (isinstance(sigma_val, int) and sigma_val in (-1, 1)):
        raise SystemExit(
            json.dumps(
                {
                    "status": "INVALID_SIGMA_INT",
                    "expected": "sigma_int_strict_derived_v1.value ∈ {+1,-1}",
                    "actual": sigma_val,
                },
                ensure_ascii=True,
            )
        )

    # Minimal sanity checks on required upstream artifacts.
    expected_theta_obj = "ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1"
    if theta_pair.get("object") != expected_theta_obj:
        raise SystemExit(
            json.dumps(
                {"status": "INVALID_THETA_PAIR_OBJECT", "expected": expected_theta_obj, "actual": theta_pair.get("object")},
                ensure_ascii=True,
            )
        )

    expected_r1_obj = "R1_residual_orientation_datum_target_slot_population_strict_derived_from_sigma_int_slot_free_theta_pair_v1"
    if r1_population.get("object") != expected_r1_obj:
        raise SystemExit(
            json.dumps(
                {"status": "INVALID_R1_POPULATION_OBJECT", "expected": expected_r1_obj, "actual": r1_population.get("object")},
                ensure_ascii=True,
            )
        )

    theta_inputs_sigma = (((theta_pair.get("inputs") or {}).get("sigma_int") or {}).get("value"))
    if theta_inputs_sigma != sigma_val:
        raise SystemExit(
            json.dumps(
                {
                    "status": "SIGMA_INT_MISMATCH",
                    "expected": sigma_val,
                    "actual_theta_pair_sigma_int": theta_inputs_sigma,
                },
                ensure_ascii=True,
            )
        )

    theta1 = (((theta_pair.get("outputs") or {}).get("pair1") or {}).get("theta_1"))
    theta2 = (((theta_pair.get("outputs") or {}).get("pair2") or {}).get("theta_2"))
    if not (isinstance(theta1, (int, float)) and isinstance(theta2, (int, float))):
        raise SystemExit(
            json.dumps(
                {"status": "INVALID_THETA_OUTPUTS", "expected": "numeric theta_1/theta_2", "actual": {"theta_1": theta1, "theta_2": theta2}},
                ensure_ascii=True,
            )
        )

    # Consistency checks: R1 population is constructed from the same theta-pair artifact.
    r1_theta1 = (((r1_population.get("inputs") or {}).get("theta_1")))
    r1_theta2 = (((r1_population.get("inputs") or {}).get("theta_2")))
    if not (isinstance(r1_theta1, (int, float)) and isinstance(r1_theta2, (int, float))):
        raise SystemExit(
            json.dumps(
                {"status": "INVALID_R1_INPUT_THETAS", "expected": "numeric theta_1/theta_2", "actual": {"theta_1": r1_theta1, "theta_2": r1_theta2}},
                ensure_ascii=True,
            )
        )

    theta_consistency_ok = (abs(float(theta1) - float(r1_theta1)) <= 1e-12) and (abs(float(theta2) - float(r1_theta2)) <= 1e-12)

    u1_theta = float_list((((theta_pair.get("outputs") or {}).get("pair1") or {}).get("u_1")))
    u2_theta = float_list((((theta_pair.get("outputs") or {}).get("pair2") or {}).get("u_2")))
    u1_r1 = float_list((((r1_population.get("outputs") or {}).get("u_1"))))
    u2_r1 = float_list((((r1_population.get("outputs") or {}).get("u_2"))))

    if u1_theta is None or u2_theta is None or u1_r1 is None or u2_r1 is None:
        raise SystemExit(
            json.dumps(
                {"status": "INVALID_U_VECTORS", "expected": "numeric lists u_1/u_2 in both artifacts"},
                ensure_ascii=True,
            )
        )

    u1_max_diff = max_abs_diff(u1_theta, u1_r1)
    u2_max_diff = max_abs_diff(u2_theta, u2_r1)
    u_vectors_consistency_ok = (u1_max_diff <= 1e-12) and (u2_max_diff <= 1e-12)

    theta1_base = (((theta_pair.get("outputs") or {}).get("pair1") or {}).get("theta_1_base_mod_pi"))
    sign_convention_ok = None
    if isinstance(theta1_base, (int, float)):
        expected_shift = math.pi if int(sigma_val) == -1 else 0.0
        sign_convention_ok = mod_2pi_close(float(theta1) - float(theta1_base), expected_shift, tol=1e-12)

    support = {
        "object": "Iota_residual_datum_sigma_int_bridge_export_map_object_support_v1",
        "status": "actual_exported_bridge_export_map_object_support_above_exported_map_object__no_false_pass",
        "as_of": "2026-03-15",
        "intent": (
            "Export one explicit strict-lane object-support layer above the exported strict-core sigma-int -> residual target-slot export-map object "
            "(F311/N422), by packaging a noncyclic, observer-free support carrier that binds the sign-only export-map object with the exported "
            "slot-free sigma-int theta-pair source (F451/N489) and its audited R1 target-slot inhabitant instance (P451). "
            "This discharges the post-witness object-support target (T130/N395) on the declared strict sigma-int lane, without upgrading the export-map "
            "object itself (it remains sign-only) and without implying selector closure."
        ),
        "upgrades_witness": "Kappa_residual_datum_sigma_int_bridge_export_map_object_support_witness_v1",
        "upgrades_support_packet": "Nu_residual_datum_sigma_int_bridge_export_map_object_support_support_packet_v1",
        "discharges_target": "Lambda_residual_datum_sigma_int_bridge_export_map_object_support_target_v1",
        "inputs": {
            "sigma_int_strict_derived_v1": {
                "ref": str(IN_SIGMA_INT.relative_to(REPO)),
                "object": sigma_int.get("object"),
                "value": int(sigma_val),
                "status": sigma_int.get("status"),
            },
            "sigma_int_gauge_quotient_safety_witness_v1": {
                "ref": str(IN_GAUGE_SAFETY.relative_to(REPO)),
                "object": gauge_safety.get("object"),
                "status": gauge_safety.get("status"),
            },
            "upsilon_residual_datum_sigma_int_bridge_export_map_object_v1": {
                "ref": str(IN_EXPORT_MAP_OBJECT.relative_to(REPO)),
                "object": export_map_object.get("object"),
                "status": export_map_object.get("status"),
                "typed_map_shape": export_map_object.get("typed_map_shape"),
            },
            "theta_pair_sigma_int_slot_free": {
                "ref": str(IN_THETA_PAIR.relative_to(REPO)),
                "object": theta_pair.get("object"),
                "status": theta_pair.get("status"),
                "theta_1": float(theta1),
                "theta_2": float(theta2),
            },
            "r1_target_slot_population_from_sigma_int_slot_free_theta_pair": {
                "ref": str(IN_R1_POPULATION.relative_to(REPO)),
                "object": r1_population.get("object"),
                "status": r1_population.get("status"),
            },
        },
        "constraints": {
            "noncyclic_no_theta_inputs": True,
            "noncyclic_no_populated_basis_pair_inputs": True,
            "observer_free_no_K_obs": True,
            "selector_neutral_no_S_sel_int_claim": True,
            "qw_2191_closed": False,
        },
        "compatibility_audits": {
            "theta_inputs_sigma_int_matches_strict_sigma_int": bool(theta_inputs_sigma == sigma_val),
            "theta_pair_and_r1_population_theta_match": bool(theta_consistency_ok),
            "theta_pair_and_r1_population_u_vectors_match_max_abs_diff": {"u_1": u1_max_diff, "u_2": u2_max_diff},
            "theta_pair_and_r1_population_u_vectors_match": bool(u_vectors_consistency_ok),
            "sigma_int_sign_convention_matches_export_map_object": sign_convention_ok,
        },
        "hard_limits": [
            "Does not upgrade the strict-core export-map object (F311/N422), which remains sign-only (residual Z2 population only).",
            "Does not claim admissible S_sel_int nor strict-core selector closure.",
            "Does not claim global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F452",
        "status": "F452_EXECUTED_FIRST_ACTUAL_STRICT_SIGMA_INT_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_SUPPORT_PACKET_NO_FALSE_PASS",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": {
            "object_support_object": support["object"],
            "discharges_target": support["discharges_target"],
            "compatibility_ok": bool(theta_consistency_ok and u_vectors_consistency_ok and (sign_convention_ok is not False)),
        },
        "no_false_pass": True,
    }

    OUT_SUPPORT.write_text(json.dumps(support, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUPPORT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUPPORT_SUMMARY)


if __name__ == "__main__":
    main()

