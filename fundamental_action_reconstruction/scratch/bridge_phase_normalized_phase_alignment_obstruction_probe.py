#!/usr/bin/env python3
"""Scratch probe: affine phase-alignment obstruction for the horizon bridge.

The compression-horizon packet found a useful candidate: the inverse branch maps
the strict carrier zero to the legacy-side horizon x_h=4/3.  A tempting stronger
claim would be that the ontological bridge is just affine phase alignment,

    omega_legacy*x + phi_legacy = omega_strict*d + phi_strict.

This probe checks that tempting claim and blocks it.  Exact phase alignment does
map the strict zero to x_h, but it sends d=0 to a negative legacy coordinate.
The affine map that preserves the origin misses the horizon, and the affine map
that preserves origin+horizon has the wrong phase slope.  Therefore any viable
ontological bridge in this lane must be nonlinear/compression-flow-like, not a
single affine phase reparameterization.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
HORIZON = HERE / "bridge_phase_normalized_compression_horizon_report.json"
OUT_JSON = HERE / "bridge_phase_normalized_phase_alignment_obstruction_report.json"
OUT_MD = HERE / "bridge_phase_normalized_phase_alignment_obstruction_report.md"

LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 9.0 / 5.0}
Z12_D = np.arange(0, 12, dtype=float)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def legacy_norm(x: np.ndarray) -> np.ndarray:
    return np.cos(LEGACY["omega"] * x + LEGACY["phi"]) / math.cos(LEGACY["phi"]) / (1.0 + LEGACY["beta_tors"] * x)


def strict_norm(d: np.ndarray) -> np.ndarray:
    return np.cos(STRICT["omega"] * d + STRICT["phi"]) / math.cos(STRICT["phi"]) / (1.0 + d ** STRICT["eta"])


def affine_maps(horizon: dict[str, Any]) -> dict[str, Any]:
    x_h = float(horizon["x_horizon"])
    d_zero = float(horizon["strict_zero_d"])
    phase_slope = STRICT["omega"] / LEGACY["omega"]
    phase_intercept = (STRICT["phi"] - LEGACY["phi"]) / LEGACY["omega"]
    origin_horizon_slope = x_h / d_zero
    origin_phase_slope = phase_slope
    return {
        "strict_carrier_zero": d_zero,
        "legacy_horizon": x_h,
        "exact_phase_alignment": {
            "map": "x_phase(d)=(omega_s*d+phi_s-phi_l)/omega_l",
            "slope": phase_slope,
            "intercept": phase_intercept,
            "x_at_d0": phase_slope * d_zero + phase_intercept,
            "x_at_0": phase_intercept,
            "horizon_residual_at_d0": phase_slope * d_zero + phase_intercept - x_h,
            "origin_admissible": phase_intercept >= 0.0,
        },
        "origin_preserving_phase_slope": {
            "map": "x_origin_phase(d)=(omega_s/omega_l)*d",
            "slope": origin_phase_slope,
            "x_at_d0": origin_phase_slope * d_zero,
            "horizon_residual_at_d0": origin_phase_slope * d_zero - x_h,
        },
        "origin_and_horizon_affine": {
            "map": "x_origin_horizon(d)=(x_h/d0)*d",
            "slope": origin_horizon_slope,
            "phase_slope_residual": LEGACY["omega"] * origin_horizon_slope - STRICT["omega"],
            "relative_slope_vs_phase_alignment": origin_horizon_slope / phase_slope,
        },
    }


def map_values(kind: str, maps: dict[str, Any], d: np.ndarray) -> np.ndarray:
    if kind == "exact_phase_alignment":
        row = maps["exact_phase_alignment"]
        return row["slope"] * d + row["intercept"]
    if kind == "origin_preserving_phase_slope":
        return maps[kind]["slope"] * d
    if kind == "origin_and_horizon_affine":
        return maps[kind]["slope"] * d
    raise ValueError(kind)


def stress_rows(maps: dict[str, Any]) -> dict[str, Any]:
    target = strict_norm(Z12_D)
    out: dict[str, Any] = {}
    for kind in ["exact_phase_alignment", "origin_preserving_phase_slope", "origin_and_horizon_affine"]:
        x = map_values(kind, maps, Z12_D)
        residual = legacy_norm(x) - target
        out[kind] = {
            "min_x_on_z12": float(np.min(x)),
            "max_x_on_z12": float(np.max(x)),
            "nonnegative_on_z12": bool(np.all(x >= 0.0)),
            "max_abs_output_residual_on_z12": float(np.max(np.abs(residual))),
            "rows": [
                {
                    "d": float(d),
                    "x": float(xi),
                    "legacy_norm_at_x": float(legacy_norm(np.array([xi]))[0]),
                    "strict_norm": float(si),
                    "residual": float(ri),
                }
                for d, xi, si, ri in zip(Z12_D, x, target, residual)
            ],
        }
    return out


def main() -> None:
    horizon_report = load_json(HORIZON)
    horizon = horizon_report["horizon"]
    maps = affine_maps(horizon)
    stress = stress_rows(maps)
    obstruction = {
        "exact_phase_alignment_maps_horizon": abs(maps["exact_phase_alignment"]["horizon_residual_at_d0"]) < 1e-12,
        "exact_phase_alignment_fails_origin_nonnegativity": not maps["exact_phase_alignment"]["origin_admissible"],
        "origin_phase_map_misses_horizon": abs(maps["origin_preserving_phase_slope"]["horizon_residual_at_d0"]) > 1e-3,
        "origin_horizon_map_has_wrong_phase_slope": abs(maps["origin_and_horizon_affine"]["phase_slope_residual"]) > 1e-3,
    }
    obstruction["all_simple_affine_routes_blocked"] = all(obstruction.values())

    report = {
        "status": "OPEN_AFFINE_PHASE_ALIGNMENT_OBSTRUCTION_NO_BRIDGE_THEOREM",
        "result_kind": "SCRATCH_PHASE_NORMALIZED_PHASE_ALIGNMENT_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        "source_reports": {"compression_horizon": str(HORIZON.relative_to(HERE.parents[1]))},
        "affine_maps": maps,
        "z12_stress": stress,
        "restricted_obstruction": obstruction,
        "candidate_ontological_reading": {
            "name": "nonlinear_compression_flow_required",
            "supported_by_this_probe": bool(obstruction["all_simple_affine_routes_blocked"]),
            "content": "The horizon bridge cannot be just affine carrier phase alignment; it requires a nonlinear compression law if it is physical.",
            "not_yet_proven_because": "The nonlinear law has not been derived from strict-side information geometry.",
        },
        "upstream_replay": {
            "horizon_candidate_supported": horizon_report["candidate_ontological_reading"]["supported_by_this_probe"],
            "horizon_x": horizon["x_horizon"],
            "strict_zero_d": horizon["strict_zero_d"],
        },
        "honest_interpretation": [
            "Strict and legacy carrier zeros can be aligned affinely, but that affine map sends d=0 to a negative legacy coordinate.",
            "The origin-preserving affine alternatives either miss the finite horizon or use the wrong phase slope.",
            "Thus the ontological bridge candidate, if real, must be a nonlinear compression flow rather than a simple phase reparameterization.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No affine phase map is accepted as a bridge.",
            "No strict-side derivation of the nonlinear compression flow is exported.",
            "No legacy physical-role transfer is licensed.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Use the affine obstruction as a constraint on future ontology: derive a nonlinear compression-flow law from strict-side geometry or classify the horizon bridge as output-matching artifact.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch phase-normalized affine phase-alignment obstruction\n\n"
        "Status: simple affine phase bridge blocked; no bridge theorem.\n\n"
        f"- Exact phase alignment maps strict zero to horizon with residual `{maps['exact_phase_alignment']['horizon_residual_at_d0']:.3e}` but has `x(0)={maps['exact_phase_alignment']['x_at_0']:.12f}`.\n"
        f"- Origin-preserving phase map misses horizon by `{maps['origin_preserving_phase_slope']['horizon_residual_at_d0']:.12f}`.\n"
        f"- Origin+horizon affine map phase-slope residual `{maps['origin_and_horizon_affine']['phase_slope_residual']:.12f}`; all simple affine routes blocked `{obstruction['all_simple_affine_routes_blocked']}`.\n"
        "- Honest read: any real ontological bridge must be nonlinear compression-flow-like, not carrier-phase affine alignment.\n"
        "- No false pass: no kernel identity, no affine bridge, no strict nonlinear-flow derivation, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
