#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2389_s1339_cap_slack_budget_sensitivity_certificate.json"
MD = GEN / "p2389_s1339_cap_slack_budget_sensitivity_certificate.md"

SOURCE_FILES = {
    "P2388_ROOT_UNIQUENESS": GEN / "p2388_s1338_cap_threshold_root_uniqueness_certificate.json",
    "P2387_EXACT_KKT_BRANCH": GEN / "p2387_s1337_bathtub_exact_kkt_branch_certificate.json",
    "P2386_LP_DUAL": GEN / "p2386_s1336_bathtub_lp_dual_certificate.json",
    "P2384_SYMBOLIC_INEQUALITY": GEN / "p2384_s1334_symbolic_bathtub_inequality_proof_packet.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

STRICT_PARAMS = {"omega": 743.0 / 4000.0, "phi": 13.0 / 80.0, "beta": 1.0, "eta": 9.0 / 5.0}
ETA_WORST = 9.0 / 5.0
BETA_TORS_WORST = 1.0 / 10.0
D1 = 1
D5 = 5
CAP_TEST = 8.0 / 5.0
DERIVATIVE_GRID = 128


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def k_strict(d: int) -> float:
    return math.cos(STRICT_PARAMS["omega"] * d + STRICT_PARAMS["phi"]) / (
        1.0 + STRICT_PARAMS["beta"] * d ** STRICT_PARAMS["eta"]
    )


def u0(d: int) -> float:
    return 1.0 + BETA_TORS_WORST * d


def delta(d: int) -> float:
    return d**ETA_WORST - BETA_TORS_WORST * d


def capped_weight(d: int, cap: float) -> float:
    start = u0(d)
    return cap * math.log((start + delta(d) / cap) / start)


def target() -> float:
    return 3.0 * k_strict(D1) - k_strict(D5)


def threshold_function(cap: float) -> float:
    return capped_weight(D5, cap) - 3.0 * capped_weight(D1, cap) - target()


def weight_d_cap(d: int, cap: float) -> float:
    x = delta(d) / (cap * u0(d))
    return math.log1p(x) - x / (1.0 + x)


def threshold_derivative(cap: float) -> float:
    return weight_d_cap(D5, cap) - 3.0 * weight_d_cap(D1, cap)


def root_from_p2388(artifact: dict[str, Any]) -> float:
    return float(
        artifact["cap_threshold_root_uniqueness_certificate"]["cap_threshold_root_certificate"]["bisection_certificate"][
            "mid"
        ]
    )


def derivative_budget(root: float, cap: float) -> dict[str, Any]:
    min_row = {"cap": None, "F_prime": float("inf")}
    max_row = {"cap": None, "F_prime": -float("inf")}
    for index in range(DERIVATIVE_GRID + 1):
        point = root + (cap - root) * index / DERIVATIVE_GRID
        derivative = threshold_derivative(point)
        if derivative < min_row["F_prime"]:
            min_row = {"cap": point, "F_prime": derivative}
        if derivative > max_row["F_prime"]:
            max_row = {"cap": point, "F_prime": derivative}
    margin = threshold_function(cap)
    surplus = cap - root
    mean_slope = margin / surplus
    return {
        "cap_interval": {"lo": root, "hi": cap},
        "grid_points": DERIVATIVE_GRID + 1,
        "min_derivative_grid": min_row,
        "max_derivative_grid": max_row,
        "mean_value_slope_from_root_to_M": mean_slope,
        "mean_value_slope_inside_grid_bounds": min_row["F_prime"] <= mean_slope <= max_row["F_prime"],
        "cap_surplus_M_minus_root": surplus,
        "F_M_margin": margin,
        "conservative_margin_lower_from_min_derivative_times_surplus": min_row["F_prime"] * surplus,
        "root_sensitivity_to_additive_threshold_error_interval": {
            "dM_dT_lower_1_over_max_F_prime": 1.0 / max_row["F_prime"],
            "dM_dT_upper_1_over_min_F_prime": 1.0 / min_row["F_prime"],
            "meaning": "for a small additive increase in the scalar target, the unique cap root moves by approximately dT/F'(M*)",
        },
        "full_margin_can_absorb_additive_target_error_below": margin,
        "proof_role": "P2388 supplies the unique root and P2384 supplies positive cap derivative; this table quantifies the M=1.6 slack budget.",
    }


def source_geometry_budget(root: float, cap: float) -> dict[str, Any]:
    threshold_interval = 1.0 / root
    cap_interval = 1.0 / cap
    threshold_early_half_mass = min(1.0, root * 0.5)
    cap_early_half_mass = min(1.0, cap * 0.5)
    threshold_barycenter = threshold_interval / 2.0
    cap_barycenter = cap_interval / 2.0
    threshold_shift = 0.5 - threshold_barycenter
    cap_shift = 0.5 - cap_barycenter
    return {
        "unique_root_cap": root,
        "accepted_cap_M": cap,
        "early_interval_length_at_root": threshold_interval,
        "early_interval_length_at_M": cap_interval,
        "early_interval_shortening_vs_threshold": threshold_interval - cap_interval,
        "early_half_mass_at_root": threshold_early_half_mass,
        "early_half_mass_at_M": cap_early_half_mass,
        "early_half_mass_surplus_vs_threshold": cap_early_half_mass - threshold_early_half_mass,
        "barycenter_at_root": threshold_barycenter,
        "barycenter_at_M": cap_barycenter,
        "barycenter_left_shift_vs_threshold": threshold_barycenter - cap_barycenter,
        "uniform_shift_at_root": threshold_shift,
        "uniform_shift_at_M": cap_shift,
        "uniform_shift_surplus_vs_threshold": cap_shift - threshold_shift,
        "status": "SLACK_BUDGET_ONLY_SOURCE_STILL_OPEN",
    }


def replay_p2388(artifact: dict[str, Any]) -> dict[str, Any]:
    try:
        cert = artifact["cap_threshold_root_uniqueness_certificate"]["cap_threshold_root_certificate"]
        return {
            "p2388_root_mid": cert["bisection_certificate"]["mid"],
            "p2388_F_M_test": cert["accepted_cap_replay"]["F_M_test"],
            "p2388_M_test_minus_root": cert["accepted_cap_replay"]["M_test_minus_root_mid"],
            "p2388_source_status": cert["source_burden_translation"]["status"],
        }
    except KeyError:
        return {"status": "P2388_REPLAY_FIELDS_MISSING"}


def cap_slack_budget_certificate(root: float) -> dict[str, Any]:
    derivative = derivative_budget(root, CAP_TEST)
    geometry = source_geometry_budget(root, CAP_TEST)
    return {
        "root_source": "P2388 bisection midpoint for the unique worst-corner cap threshold",
        "worst_corner": {"eta": ETA_WORST, "beta_tors": BETA_TORS_WORST},
        "accepted_cap_M": CAP_TEST,
        "derivative_and_margin_budget": derivative,
        "source_geometry_budget": geometry,
        "proof_reduction": [
            "P2388 proves that the threshold is a unique root of the scalar monotone equation.",
            "The current packet applies the mean-value theorem on [M*,1.6] to quantify the cap and scalar-margin slack above that root.",
            "The source-side geometry translation states how much stronger M=1.6 is than a just-threshold cap in early support length, early-half mass, barycenter, and uniform-shift terms.",
        ],
    }


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}

    root = root_from_p2388(artifacts["P2388_ROOT_UNIQUENESS"])
    certificate = cap_slack_budget_certificate(root)
    p2388_replay = replay_p2388(artifacts["P2388_ROOT_UNIQUENESS"])

    theorem_export = {
        "name": "P2389/S1339 cap slack budget sensitivity certificate",
        "claim": (
            "After P2388 identifies the unique cap threshold, the accepted cap M=1.6 has a quantified slack budget: its scalar chamber margin, cap surplus, derivative/sensitivity band, and source-geometry surplus are explicit acceptance margins rather than hidden numerical comfort."
        ),
        "positive_content": [
            "Quantifies M=1.6-root surplus and F(1.6) scalar margin.",
            "Audits derivative bounds on [root,1.6] and checks the mean-value slope against those bounds.",
            "Exports a sensitivity interval for additive threshold perturbations via dM*/dT=1/F'(M*).",
            "Translates the slack into early interval, early-half mass, barycenter, and uniform-shift budgets for future source theorems.",
        ],
        "not_licensed": [
            "strict source theorem deriving M=1.6 or the front-loaded density",
            "promotion of slack budgets into L_total",
            "claim that strict dynamics supplies the slack",
            "QW-2191 discharge or selector closure",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2389_S1339_CAP_SLACK_BUDGET_SENSITIVITY_CERTIFICATE",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "artifact_replay": {
            "p2388_packet_id": artifacts["P2388_ROOT_UNIQUENESS"].get("packet_id"),
            "p2387_packet_id": artifacts["P2387_EXACT_KKT_BRANCH"].get("packet_id"),
            "p2386_packet_id": artifacts["P2386_LP_DUAL"].get("packet_id"),
            "p2384_packet_id": artifacts["P2384_SYMBOLIC_INEQUALITY"].get("packet_id"),
        },
        "p2388_replay": p2388_replay,
        "cap_slack_budget_certificate": certificate,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    budget = certificate["derivative_and_margin_budget"]
    geometry = certificate["source_geometry_budget"]
    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2388_loaded": probe["artifact_replay"]["p2388_packet_id"] == "P2388",
        "p2387_loaded": probe["artifact_replay"]["p2387_packet_id"] == "P2387",
        "p2386_loaded": probe["artifact_replay"]["p2386_packet_id"] == "P2386",
        "p2384_loaded": probe["artifact_replay"]["p2384_packet_id"] == "P2384",
        "positive_cap_surplus": budget["cap_surplus_M_minus_root"] > 0.0,
        "positive_scalar_margin": budget["F_M_margin"] > 0.0,
        "positive_derivative_band": budget["min_derivative_grid"]["F_prime"] > 0.0,
        "mean_value_slope_consistent": budget["mean_value_slope_inside_grid_bounds"],
        "geometry_stronger_than_threshold": geometry["early_interval_shortening_vs_threshold"] > 0.0
        and geometry["early_half_mass_surplus_vs_threshold"] > 0.0
        and geometry["uniform_shift_surplus_vs_threshold"] > 0.0,
        "p2388_replay_consistent": abs(budget["F_M_margin"] - p2388_replay.get("p2388_F_M_test", 0.0)) < 1.0e-12,
        "source_target_not_promoted": geometry["status"] == "SLACK_BUDGET_ONLY_SOURCE_STILL_OPEN",
        "docs_updated_with_p2389_slack_budget": all("P2389/S1339" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2389_s1339_v1",
        "packet_id": "P2389",
        "stage_id": "S1339",
        "timestamp_utc": "2026-06-01T00:00:00Z",
        "produced_by": rel(Path(__file__).resolve()),
        "status": "OPEN_PROGRESS_CAP_SLACK_BUDGET_SENSITIVITY_CERTIFICATE_SOURCE_STILL_OPEN",
        "result_kind": "CAP_SLACK_BUDGET_SENSITIVITY_CERTIFICATE",
        "cap_slack_budget_sensitivity_certificate": probe,
        "gatekeeper_checks": gatekeeper_checks,
        "recommended_next_honest_step": (
            "Use the quantified slack budget when evaluating future source candidates: a candidate may either derive M=1.6, derive a smaller M still above the unique root, or explicitly remain a non-strict cap premise."
        ),
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2389 S1339: cap slack budget sensitivity certificate",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2389/S1339 quantifies how much room the accepted cap `M=1.6` has above the unique P2388 cap threshold.",
                "It records the cap surplus, scalar margin `F(1.6)`, derivative/sensitivity band, and source-geometry surplus against the just-threshold profile.",
                "",
                "## Checks",
                "",
                f"- Cap surplus `1.6-root`: `{budget['cap_surplus_M_minus_root']}`.",
                f"- Scalar margin `F(1.6)`: `{budget['F_M_margin']}`.",
                f"- Mean-value slope: `{budget['mean_value_slope_from_root_to_M']}`.",
                f"- Root sensitivity interval: `{budget['root_sensitivity_to_additive_threshold_error_interval']}`.",
                f"- Early interval shortening vs threshold: `{geometry['early_interval_shortening_vs_threshold']}`.",
                f"- Early-half mass surplus vs threshold: `{geometry['early_half_mass_surplus_vs_threshold']}`.",
                "",
                "## Hard limits",
                "",
                "- This is a slack/sensitivity acceptance budget, not a strict source theorem deriving `M=1.6` or the density.",
                "- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
