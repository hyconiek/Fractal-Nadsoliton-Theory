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

IN_THETA_PAIR = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"

OUT_ALPHA12 = (
    GENERATED
    / "alpha12_pair1_pair2_transition_angle_strict_derived_from_sigma_int_slot_free_theta_pair_v1.json"
)
OUT_ALPHA12_SUMMARY = (
    GENERATED
    / "alpha12_pair1_pair2_transition_angle_strict_derived_from_sigma_int_slot_free_theta_pair_v1_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def angle_mod(x: float, modulus: float) -> float:
    return float(x % modulus)


def signed_angle_from_mod_2pi(x_mod_2pi: float) -> float:
    # Map [0,2π) -> (-π,π]
    if x_mod_2pi > math.pi:
        return float(x_mod_2pi - 2.0 * math.pi)
    return float(x_mod_2pi)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_THETA_PAIR.exists():
        raise SystemExit(
            json.dumps(
                {
                    "status": "NOT_COMPUTABLE_MISSING_INPUT",
                    "missing": str(IN_THETA_PAIR.relative_to(REPO)),
                    "expected": "F451 theta-pair export object",
                },
                ensure_ascii=True,
            )
        )

    theta_pair = load_json(IN_THETA_PAIR)
    try:
        theta_1 = float(((theta_pair.get("outputs") or {}).get("pair1") or {}).get("theta_1"))
        theta_2 = float(((theta_pair.get("outputs") or {}).get("pair2") or {}).get("theta_2"))
        sigma_int = int(((theta_pair.get("inputs") or {}).get("sigma_int") or {}).get("value"))
    except Exception as exc:
        raise SystemExit(
            json.dumps(
                {
                    "status": "INVALID_INPUT_SHAPE",
                    "expected": "F451.outputs.pair1.theta_1, F451.outputs.pair2.theta_2, F451.inputs.sigma_int.value",
                    "error": str(exc),
                },
                ensure_ascii=True,
            )
        )

    alpha_12_raw = float(theta_2 - theta_1)
    alpha_12_mod_2pi = angle_mod(alpha_12_raw, 2.0 * math.pi)
    alpha_12_mod_pi = angle_mod(alpha_12_raw, math.pi)
    alpha_12_signed = signed_angle_from_mod_2pi(alpha_12_mod_2pi)

    c = math.cos(alpha_12_mod_2pi)
    s = math.sin(alpha_12_mod_2pi)
    g12 = [[c, -s], [s, c]]

    artifact = {
        "object": "Alpha12_pair1_pair2_transition_angle_strict_derived_from_sigma_int_slot_free_theta_pair_v1",
        "status": "actual_exported_strict_derived_transition_angle_alpha12__pair1_pair2__no_false_pass",
        "as_of": "2026-03-15",
        "intent": (
            "Export an explicit strict-derived pair1/pair2 transition angle alpha_12 := (theta_2 - theta_1) mod 2π "
            "derived only from the already exported strict-core sigma-int slot-free theta-pair supply (F451/N489). "
            "This is lane-scoped and does not imply any global selector transition/gluing object."
        ),
        "inputs": {
            "theta_pair_source": {
                "object": "ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1",
                "ref": str(IN_THETA_PAIR.relative_to(REPO)),
                "sigma_int": sigma_int,
                "theta_1": theta_1,
                "theta_2": theta_2,
            }
        },
        "outputs": {
            "alpha_12_raw": alpha_12_raw,
            "alpha_12_mod_2pi": alpha_12_mod_2pi,
            "alpha_12_signed": alpha_12_signed,
            "alpha_12_mod_pi": alpha_12_mod_pi,
            "G12_orbit_coordinate_rotation_mod_2pi": g12,
        },
        "hard_limits": [
            "Does not claim any global discharge of QW-2191.",
            "Does not claim admissible S_sel_int nor strict-core selector closure.",
            "Does not claim a global selector transition/gluing object beyond the declared sigma-int lane.",
            "Does not assign a theorem-level physical interpretation to alpha_12.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F457",
        "status": "F457_EXECUTED_CURRENT_STRICT_SIGMA_INT_THETA_PAIR_TO_ALPHA12_TRANSITION_ANGLE_EXPORT_PACKET_NO_FALSE_PASS",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "theta_pair_ref": str(IN_THETA_PAIR.relative_to(REPO)),
            "sigma_int": sigma_int,
        },
        "outputs": {
            "theta_1": theta_1,
            "theta_2": theta_2,
            "alpha_12_mod_2pi": alpha_12_mod_2pi,
            "alpha_12_mod_pi": alpha_12_mod_pi,
        },
        "no_false_pass": True,
    }

    OUT_ALPHA12.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_ALPHA12_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_ALPHA12_SUMMARY)


if __name__ == "__main__":
    main()

