#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from itertools import combinations, product
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2381_s1331_affine_frontload_source_obligation_theorem.json"
MD = GEN / "p2381_s1331_affine_frontload_source_obligation_theorem.md"

SOURCE_FILES = {
    "P2380_RECTANGLE_MONOTONICITY": GEN / "p2380_s1330_front_loaded_profile_rectangle_monotonicity_certificate.json",
    "P2379_FRONT_LOADED_PROFILE": GEN / "p2379_s1329_front_loaded_normalized_transport_density_existence_probe.json",
    "P2378_UNIT_NORMALIZED_INSUFFICIENCY": GEN / "p2378_s1328_unit_normalized_transport_coupling_insufficiency_theorem.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
CHAMBER_THRESHOLD = 1.0 / 3.0
ETA_INTERVAL = (9.0 / 5.0, 2.0)
BETA_TORS_INTERVAL = (0.0, 0.1)
SAMPLE_ETAS = [9.0 / 5.0, 19.0 / 10.0, 2.0]
SAMPLE_BETA_TORS = [0.0, 0.01, 0.1]
STRICT_BETA = 1.0
STRICT_PARAMS = {"omega": 743.0 / 4000.0, "phi": 13.0 / 80.0, "beta": STRICT_BETA, "eta": 9.0 / 5.0}
LAMBDA_TEST = 4.0 / 5.0
EPSILON = 1e-9


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


def compression_log_weight(d: int, eta: float, beta_tors: float) -> float:
    return math.log((1.0 + STRICT_BETA * d**eta) / (1.0 + beta_tors * d))


def affine_frontload_moment(d: int, eta: float, beta_tors: float) -> float:
    start = 1.0 + beta_tors * d
    delta = d**eta - beta_tors * d
    c = compression_log_weight(d, eta, beta_tors)
    return c * (1.0 + 2.0 * start / delta) - 2.0


def affine_weighted_transport(d: int, eta: float, beta_tors: float, lam: float) -> float:
    return compression_log_weight(d, eta, beta_tors) + lam * affine_frontload_moment(d, eta, beta_tors)


def lambda_threshold(eta: float, beta_tors: float) -> dict[str, float]:
    k1 = k_strict(1)
    k5 = k_strict(5)
    c1 = compression_log_weight(1, eta, beta_tors)
    c5 = compression_log_weight(5, eta, beta_tors)
    b1 = affine_frontload_moment(1, eta, beta_tors)
    b5 = affine_frontload_moment(5, eta, beta_tors)
    numerator = 3.0 * (k1 + c1) - (k5 + c5)
    denominator = b5 - 3.0 * b1
    return {
        "K1": k1,
        "K5": k5,
        "C1": c1,
        "C5": c5,
        "B1": b1,
        "B5": b5,
        "lambda_numerator": numerator,
        "lambda_denominator_B5_minus_3B1": denominator,
        "lambda_threshold_gt": numerator / denominator,
    }


def affine_profile_obligations(lam: float) -> dict[str, float]:
    return {
        "lambda": lam,
        "rho_s0": 1.0 + lam,
        "rho_s1": 1.0 - lam,
        "endpoint_contrast_rho0_over_rho1": (1.0 + lam) / (1.0 - lam),
        "early_half_mass_int_0_to_1_2": 0.5 + lam / 4.0,
        "late_half_mass_int_1_2_to_1": 0.5 - lam / 4.0,
        "transport_time_barycenter_int_s_rho": 0.5 - lam / 6.0,
        "barycenter_shift_from_uniform": lam / 6.0,
        "l1_distance_from_uniform": lam / 2.0,
        "l2_squared_distance_from_uniform": lam * lam / 3.0,
        "profile_variance_int_rho_minus_1_sq": lam * lam / 3.0,
    }


def internal_edges(support: tuple[int, ...], step: int) -> int:
    support_set = set(support)
    edges = set()
    for node in support_set:
        for neighbor in ((node + step) % Z12_NODE_COUNT, (node - step) % Z12_NODE_COUNT):
            if neighbor in support_set:
                edges.add(tuple(sorted((node, neighbor))))
    return len(edges)


def support_rows() -> list[dict[str, Any]]:
    return [
        {"support": list(support), "h1": internal_edges(support, 1), "h5": internal_edges(support, 5)}
        for support in combinations(range(Z12_NODE_COUNT), SUPPORT_SIZE)
    ]


def score_maximizers(rows: list[dict[str, Any]], a: float, b: float) -> dict[str, Any]:
    scored = [(a * row["h1"] + b * row["h5"], row) for row in rows]
    maximum = max(score for score, _ in scored)
    maximizers = [row for score, row in scored if abs(score - maximum) <= 1e-10]
    pair_distribution: dict[str, int] = {}
    for row in maximizers:
        key = f"{row['h1']},{row['h5']}"
        pair_distribution[key] = pair_distribution.get(key, 0) + 1
    return {
        "a_over_b": a / b if b > 0 else None,
        "d5_chamber": b > 0 and a >= 0 and a / b < CHAMBER_THRESHOLD,
        "maximum_score": maximum,
        "maximizer_count": len(maximizers),
        "maximizer_h1_h5_pair_distribution": dict(sorted(pair_distribution.items())),
    }


def score_for_lambda(eta: float, beta_tors: float, lam: float) -> dict[str, Any]:
    threshold = lambda_threshold(eta, beta_tors)
    w1 = affine_weighted_transport(1, eta, beta_tors, lam)
    w5 = affine_weighted_transport(5, eta, beta_tors, lam)
    return score_maximizers(support_rows(), threshold["K1"] + w1, threshold["K5"] + w5)


def source_obligation_certificate() -> dict[str, Any]:
    corner_eta = ETA_INTERVAL[0]
    corner_beta_tors = BETA_TORS_INTERVAL[1]
    corner = lambda_threshold(corner_eta, corner_beta_tors)
    lambda_required = corner["lambda_threshold_gt"]
    threshold_minus = lambda_required - EPSILON
    sample_rows = []
    for eta, beta_tors in product(SAMPLE_ETAS, SAMPLE_BETA_TORS):
        row = lambda_threshold(eta, beta_tors)
        sample_rows.append(
            {
                "eta": eta,
                "beta_tors": beta_tors,
                "lambda_threshold_gt": row["lambda_threshold_gt"],
                "necessary_obligations_at_local_threshold": affine_profile_obligations(row["lambda_threshold_gt"]),
                "lambda_test_margin": LAMBDA_TEST - row["lambda_threshold_gt"],
            }
        )
    return {
        "affine_profile_family": {
            "rho_lambda": "rho_lambda(s)=1+lambda*(1-2s), 0<=lambda<1",
            "normalization": "int_0^1 rho_lambda(s) ds = 1",
            "positivity": "rho_lambda(s)>=0 iff lambda<=1",
            "acceptance_condition": "lambda > T(eta,beta_tors)",
        },
        "rectangle_worst_case_source_obligation": {
            "worst_corner": {"eta": corner_eta, "beta_tors": corner_beta_tors, **corner},
            "lambda_must_exceed": lambda_required,
            "strict_inequality_note": "At the threshold the score lies on the chamber wall; below it the d5 chamber fails at the worst corner.",
            "necessary_profile_obligations_for_any_rectangle_uniform_affine_source": affine_profile_obligations(lambda_required),
            "lambda_test_obligations": affine_profile_obligations(LAMBDA_TEST),
            "below_threshold_negative_control_lambda": threshold_minus,
            "below_threshold_negative_control_score": score_for_lambda(corner_eta, corner_beta_tors, threshold_minus),
            "lambda_test_corner_score": score_for_lambda(corner_eta, corner_beta_tors, LAMBDA_TEST),
        },
        "sample_local_obligation_table": sample_rows,
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}
    certificate = source_obligation_certificate()
    obligation = certificate["rectangle_worst_case_source_obligation"]
    sample_rows = certificate["sample_local_obligation_table"]

    theorem_export = {
        "theorem_name": "P2381 affine frontload source obligation theorem",
        "claim": (
            "P2380 makes lambda=0.8 sufficient for the affine normalized profile family, but P2381 records the necessary source burden. "
            "Any rectangle-uniform affine-profile source must derive lambda greater than the worst-corner threshold 0.7916644842269429, "
            "equivalently a transport barycenter at most about 0.3680559, early-half mass at least about 0.6979161, and endpoint contrast above about 8.5995."
        ),
        "positive_content": [
            "Converts the P2380 threshold into necessary affine-density asymmetry obligations.",
            "Provides a below-threshold negative control at the worst corner and a lambda=0.8 positive replay over all 792 supports.",
            "Records local source-obligation values on the P2376 3x3 grid.",
        ],
        "not_licensed": [
            "strict variational source theorem deriving the required affine asymmetry",
            "promotion of rho_lambda or weighted transport into L_total",
            "claim that arbitrary normalized profiles are solved",
            "QW-2191 discharge",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2381_S1331_AFFINE_FRONTLOAD_SOURCE_OBLIGATION_THEOREM",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "artifact_replay": {
            "p2380_packet_id": artifacts["P2380_RECTANGLE_MONOTONICITY"].get("packet_id"),
            "p2379_packet_id": artifacts["P2379_FRONT_LOADED_PROFILE"].get("packet_id"),
            "p2378_packet_id": artifacts["P2378_UNIT_NORMALIZED_INSUFFICIENCY"].get("packet_id"),
        },
        "affine_frontload_source_obligation_certificate": certificate,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2380_loaded": probe["artifact_replay"]["p2380_packet_id"] == "P2380",
        "p2379_loaded": probe["artifact_replay"]["p2379_packet_id"] == "P2379",
        "p2378_loaded": probe["artifact_replay"]["p2378_packet_id"] == "P2378",
        "worst_case_threshold_below_lambda_test": obligation["lambda_must_exceed"] < LAMBDA_TEST,
        "below_threshold_negative_control_fails_d5": not obligation["below_threshold_negative_control_score"]["d5_chamber"],
        "lambda_test_corner_selects_d5": obligation["lambda_test_corner_score"]["maximizer_h1_h5_pair_distribution"] == {"0,4": 12},
        "all_sample_thresholds_below_lambda_test": all(row["lambda_threshold_gt"] < LAMBDA_TEST for row in sample_rows),
        "docs_updated_with_p2381_source_obligation": all("P2381/S1331" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2381_s1331_v1",
        "packet_id": "P2381",
        "stage_id": "S1331",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_AFFINE_FRONTLOAD_SOURCE_OBLIGATION_QUANTIFIED_SOURCE_STILL_OPEN",
        "result_kind": "AFFINE_FRONTLOAD_SOURCE_OBLIGATION_THEOREM",
        "affine_frontload_source_obligation_theorem": probe,
        "recommended_next_honest_step": (
            "Do not treat the P2379/P2380 affine profile as sourced merely because it is sufficient. "
            "A future source theorem must derive at least the P2381 asymmetry burden, or the affine profile remains an explicit non-strict premise."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2381 S1331: affine frontload source obligation theorem",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2381/S1331 turns the P2380 sufficient affine front-loaded profile into a quantified source-obligation statement.",
                "Any rectangle-uniform affine source must derive `lambda` above the worst-corner threshold, not merely name a normalized profile.",
                "",
                "## Certificate",
                "",
                f"- Worst-corner lambda threshold: `{obligation['lambda_must_exceed']}`.",
                f"- Necessary obligations at that threshold: `{obligation['necessary_profile_obligations_for_any_rectangle_uniform_affine_source']}`.",
                f"- Lambda=0.8 obligations: `{obligation['lambda_test_obligations']}`.",
                f"- Below-threshold negative control selects: `{obligation['below_threshold_negative_control_score']['maximizer_h1_h5_pair_distribution']}`.",
                f"- Lambda=0.8 corner replay selects: `{obligation['lambda_test_corner_score']['maximizer_h1_h5_pair_distribution']}`.",
                "",
                "## Hard limits",
                "",
                "- This is a necessary source-obligation theorem for the affine family, not a derivation of that source from strict dynamics.",
                "- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
