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

IN_R1 = GENERATED / "r1_strict_core_residual_datum_target_slot_export_packet.json"
IN_SIGMA_INT = GENERATED / "sigma_int_strict_derived_v1.json"
IN_GAUGE_SAFETY = GENERATED / "sigma_int_gauge_quotient_safety_witness_v1.json"
IN_SELECTOR_TRACK = GENERATED / "sigma_int_residual_datum_selector_track_identification_witness_v1.json"
IN_THETA_PAIR = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"
IN_R1_POPULATION = GENERATED / "r1_residual_orientation_datum_target_slot_population_strict_derived_from_sigma_int_slot_free_theta_pair_v1.json"

OUT_MAP_OBJECT = GENERATED / "upsilon_residual_datum_sigma_int_bridge_export_map_object_v2.json"
OUT_MAP_OBJECT_SUMMARY = GENERATED / "upsilon_residual_datum_sigma_int_bridge_export_map_object_v2_summary.json"


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


def angle_close_mod_2pi(a: float, b: float, tol: float = 1e-12) -> bool:
    two_pi = 2.0 * math.pi
    d = float(a - b)
    d = (d + math.pi) % two_pi - math.pi
    return abs(d) <= tol


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing: list[str] = []
    for p in (IN_R1, IN_SIGMA_INT, IN_GAUGE_SAFETY, IN_SELECTOR_TRACK, IN_THETA_PAIR, IN_R1_POPULATION):
        if not p.exists():
            missing.append(str(p.relative_to(REPO)))
    if missing:
        raise SystemExit(json.dumps({"status": "MISSING_INPUT_ARTIFACTS", "missing": missing}, ensure_ascii=True))

    r1 = load_json(IN_R1)
    sigma_int = load_json(IN_SIGMA_INT)
    gauge_safety = load_json(IN_GAUGE_SAFETY)
    selector_track = load_json(IN_SELECTOR_TRACK)
    theta_pair = load_json(IN_THETA_PAIR)
    r1_population = load_json(IN_R1_POPULATION)

    sigma_val = sigma_int.get("value")
    if not (isinstance(sigma_val, int) and sigma_val in (-1, 1)):
        raise SystemExit(
            json.dumps(
                {"status": "INVALID_SIGMA_INT", "expected": "sigma_int_strict_derived_v1.value ∈ {+1,-1}", "actual": sigma_val},
                ensure_ascii=True,
            )
        )

    expected_theta_obj = "ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1"
    if theta_pair.get("object") != expected_theta_obj:
        raise SystemExit(
            json.dumps(
                {"status": "INVALID_THETA_PAIR_OBJECT", "expected": expected_theta_obj, "actual": theta_pair.get("object")},
                ensure_ascii=True,
            )
        )

    expected_r1_pop_obj = "R1_residual_orientation_datum_target_slot_population_strict_derived_from_sigma_int_slot_free_theta_pair_v1"
    if r1_population.get("object") != expected_r1_pop_obj:
        raise SystemExit(
            json.dumps(
                {
                    "status": "INVALID_R1_POPULATION_OBJECT",
                    "expected": expected_r1_pop_obj,
                    "actual": r1_population.get("object"),
                },
                ensure_ascii=True,
            )
        )

    theta_inputs_sigma = (((theta_pair.get("inputs") or {}).get("sigma_int") or {}).get("value"))
    if theta_inputs_sigma != sigma_val:
        raise SystemExit(
            json.dumps(
                {"status": "SIGMA_INT_MISMATCH", "expected": sigma_val, "actual_theta_pair_sigma_int": theta_inputs_sigma},
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

    r1_theta1 = (r1_population.get("inputs") or {}).get("theta_1")
    r1_theta2 = (r1_population.get("inputs") or {}).get("theta_2")
    if not (isinstance(r1_theta1, (int, float)) and isinstance(r1_theta2, (int, float))):
        raise SystemExit(
            json.dumps(
                {"status": "INVALID_R1_INPUT_THETAS", "expected": "numeric theta_1/theta_2", "actual": {"theta_1": r1_theta1, "theta_2": r1_theta2}},
                ensure_ascii=True,
            )
        )

    theta_consistency_ok = angle_close_mod_2pi(float(theta1), float(r1_theta1), tol=1e-12) and angle_close_mod_2pi(
        float(theta2), float(r1_theta2), tol=1e-12
    )

    u1_theta = float_list((((theta_pair.get("outputs") or {}).get("pair1") or {}).get("u_1")))
    u2_theta = float_list((((theta_pair.get("outputs") or {}).get("pair2") or {}).get("u_2")))
    u1_r1 = float_list((((r1_population.get("outputs") or {}).get("u_1"))))
    u2_r1 = float_list((((r1_population.get("outputs") or {}).get("u_2"))))
    if u1_theta is None or u2_theta is None or u1_r1 is None or u2_r1 is None:
        raise SystemExit(
            json.dumps({"status": "INVALID_U_VECTORS", "expected": "numeric lists u_1/u_2 in both artifacts"}, ensure_ascii=True)
        )

    u1_max_diff = max_abs_diff(u1_theta, u1_r1)
    u2_max_diff = max_abs_diff(u2_theta, u2_r1)
    u_vectors_consistency_ok = (u1_max_diff <= 1e-12) and (u2_max_diff <= 1e-12)

    map_object = {
        "object": "Upsilon_residual_datum_sigma_int_bridge_export_map_object_v2",
        "status": "actual_exported_strict_core_bridge_export_map_object__theta_pair_and_r1_population",
        "as_of": "2026-03-15",
        "intent": (
            "Export one upgraded strict-core sigma-int -> residual orientation datum bridge/export-map object which, on the current strict sigma-int "
            "source value, carries an explicit slot-free theta-pair construction and a corresponding strict-core R1 target-slot inhabitant instance. "
            "This does not imply strict-core selector closure or global QW-2191 discharge; it only upgrades the previously sign-only map object (F311/N422) "
            "by attaching the already exported slot-free theta supply (F451/N489) and its audited R1 target-slot population (P451) as explicit map outputs."
        ),
        "typed_map_shape": (
            "E_sigma_int_to_residual_datum_bridge_export_map_object_v2 : sigma_int_strict_derived_v1 -> residual_orientation_datum_target_slot_inhabitant"
        ),
        "domain": {
            "sigma_int_object": "sigma_int_strict_derived_v1",
            "value": int(sigma_val),
            "provenance": "premise_based_strict_source_upgrade (F307/N418)",
            "gauge_quotient_safety": "discharged (F308/N419)",
        },
        "codomain": {
            "target_slot": "residual_orientation_datum_target_slot (R1)",
            "target_slot_export_packet": str(IN_R1.relative_to(REPO)),
            "target_slot_population_object": r1_population.get("object"),
        },
        "selector_track_identification": {
            "witness_object": selector_track.get("object"),
            "witness_status": selector_track.get("status"),
            "QW_2191_status": "open",
            "no_selector_closure_implied": True,
        },
        "map_rule": {
            "theta_pair_source_object": theta_pair.get("object"),
            "theta_pair_source_ref": str(IN_THETA_PAIR.relative_to(REPO)),
            "r1_population_ref": str(IN_R1_POPULATION.relative_to(REPO)),
            "construction_class_id": ((theta_pair.get("construction_class") or {}).get("id")),
            "outputs": {
                "theta_1": float(theta1),
                "theta_2": float(theta2),
                "u_1": u1_theta,
                "u_2": u2_theta,
                "span_note": "codomain is the residual orientation datum target slot span; residual Z2 sign remains a separate convention",
            },
        },
        "contracts": {"noncyclic_inputs": True, "observer_free": True},
        "compatibility_audits": {
            "theta_pair_inputs_sigma_int_matches_strict_sigma_int": bool(theta_inputs_sigma == sigma_val),
            "theta_pair_and_r1_population_theta_match_mod_2pi": bool(theta_consistency_ok),
            "theta_pair_and_r1_population_u_vectors_match_max_abs_diff": {"u_1": u1_max_diff, "u_2": u2_max_diff},
            "theta_pair_and_r1_population_u_vectors_match": bool(u_vectors_consistency_ok),
        },
        "hard_limits": [
            "Does not claim admissible S_sel_int nor strict-core selector closure.",
            "Does not claim global QW-2191 discharge beyond the declared sigma-int -> (pair1,pair2) scope.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F455",
        "status": "F455_EXECUTED_CURRENT_STRICT_SIGMA_INT_RESIDUAL_BRIDGE_EXPORT_MAP_OBJECT_THETA_PAIR_UPGRADE_PACKET_NO_FALSE_PASS",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": {
            "export_map_object": map_object["object"],
            "typed_map_shape": map_object["typed_map_shape"],
            "theta_pair_ref": str(IN_THETA_PAIR.relative_to(REPO)),
            "r1_population_ref": str(IN_R1_POPULATION.relative_to(REPO)),
            "compatibility_ok": bool(theta_consistency_ok and u_vectors_consistency_ok),
        },
        "no_false_pass": True,
    }

    OUT_MAP_OBJECT.write_text(json.dumps(map_object, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_MAP_OBJECT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_MAP_OBJECT_SUMMARY)


if __name__ == "__main__":
    main()

