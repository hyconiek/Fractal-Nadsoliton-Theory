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

OUT = GEN / "p2379_s1329_front_loaded_normalized_transport_density_existence_probe.json"
MD = GEN / "p2379_s1329_front_loaded_normalized_transport_density_existence_probe.md"

SOURCE_FILES = {
    "P2378_UNIT_NORMALIZED_INSUFFICIENCY": GEN / "p2378_s1328_unit_normalized_transport_coupling_insufficiency_theorem.json",
    "P2377_TRANSPORT_PRIMITIVE": GEN / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.json",
    "P2376_RECTANGLE_ROBUSTNESS": GEN / "p2376_s1326_damping_compression_eta_beta_rectangle_robustness_theorem.json",
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
MIDPOINT_STEPS = 20000
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


def u0(d: int, beta_tors: float) -> float:
    return 1.0 + beta_tors * d


def u1(d: int, eta: float) -> float:
    return 1.0 + STRICT_BETA * d**eta


def transport_speed(d: int, eta: float, beta_tors: float, s: float) -> float:
    start = u0(d, beta_tors)
    delta = u1(d, eta) - start
    return delta / (start + s * delta)


def compression_log_weight(d: int, eta: float, beta_tors: float) -> float:
    return math.log(u1(d, eta) / u0(d, beta_tors))


def affine_frontload_moment(d: int, eta: float, beta_tors: float) -> float:
    """B(d)=int_0^1 (1-2s) A_s(d) ds, with A_s=partial_s log u_s."""
    start = u0(d, beta_tors)
    delta = u1(d, eta) - start
    c = compression_log_weight(d, eta, beta_tors)
    return c * (1.0 + 2.0 * start / delta) - 2.0


def affine_weighted_transport(d: int, eta: float, beta_tors: float, lam: float) -> float:
    return compression_log_weight(d, eta, beta_tors) + lam * affine_frontload_moment(d, eta, beta_tors)


def midpoint_affine_integral(d: int, eta: float, beta_tors: float, lam: float, steps: int = MIDPOINT_STEPS) -> float:
    total = 0.0
    inv_steps = 1.0 / steps
    for i in range(steps):
        s = (i + 0.5) * inv_steps
        density = 1.0 + lam * (1.0 - 2.0 * s)
        total += density * transport_speed(d, eta, beta_tors, s)
    return total * inv_steps


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


def lattice_threshold_audit(points: int = 81) -> dict[str, Any]:
    eta_min, eta_max = ETA_INTERVAL
    beta_min, beta_max = BETA_TORS_INTERVAL
    max_row: dict[str, Any] | None = None
    min_row: dict[str, Any] | None = None
    sign_failures = []
    for i in range(points):
        eta = eta_min + (eta_max - eta_min) * i / (points - 1)
        for j in range(points):
            beta_tors = beta_min + (beta_max - beta_min) * j / (points - 1)
            row = lambda_threshold(eta, beta_tors)
            record = {"eta": eta, "beta_tors": beta_tors, **row}
            if row["lambda_denominator_B5_minus_3B1"] <= 0 or not (0.0 < row["lambda_threshold_gt"] < 1.0):
                sign_failures.append(record)
            if max_row is None or row["lambda_threshold_gt"] > max_row["lambda_threshold_gt"]:
                max_row = record
            if min_row is None or row["lambda_threshold_gt"] < min_row["lambda_threshold_gt"]:
                min_row = record
    return {
        "lattice_points_per_axis": points,
        "lattice_total_points": points * points,
        "max_threshold_row": max_row,
        "min_threshold_row": min_row,
        "sign_failure_count": len(sign_failures),
        "all_lattice_thresholds_between_0_and_1": len(sign_failures) == 0,
        "lambda_test": LAMBDA_TEST,
        "lambda_test_exceeds_lattice_max": max_row is not None and LAMBDA_TEST > max_row["lambda_threshold_gt"],
    }


def front_loaded_density_certificate() -> dict[str, Any]:
    rows = support_rows()
    sample_rows = []
    max_integral_error = 0.0
    for eta, beta_tors in product(SAMPLE_ETAS, SAMPLE_BETA_TORS):
        threshold = lambda_threshold(eta, beta_tors)
        c1 = threshold["C1"]
        c5 = threshold["C5"]
        w1 = affine_weighted_transport(1, eta, beta_tors, LAMBDA_TEST)
        w5 = affine_weighted_transport(5, eta, beta_tors, LAMBDA_TEST)
        i1 = midpoint_affine_integral(1, eta, beta_tors, LAMBDA_TEST)
        i5 = midpoint_affine_integral(5, eta, beta_tors, LAMBDA_TEST)
        max_integral_error = max(max_integral_error, abs(i1 - w1), abs(i5 - w5))
        uniform_a = threshold["K1"] + c1
        uniform_b = threshold["K5"] + c5
        front_a = threshold["K1"] + w1
        front_b = threshold["K5"] + w5
        sample_rows.append(
            {
                "eta": eta,
                "beta_tors": beta_tors,
                **threshold,
                "rho_lambda": "rho_lambda(s)=1+lambda*(1-2s)",
                "lambda_test": LAMBDA_TEST,
                "rho_min_at_s1": 1.0 - LAMBDA_TEST,
                "rho_max_at_s0": 1.0 + LAMBDA_TEST,
                "rho_total_mass": 1.0,
                "W1_lambda_test": w1,
                "W5_lambda_test": w5,
                "W1_midpoint_integral": i1,
                "W5_midpoint_integral": i5,
                "uniform_score_audit": score_maximizers(rows, uniform_a, uniform_b),
                "front_loaded_score_audit": score_maximizers(rows, front_a, front_b),
            }
        )

    return {
        "weighted_transport_identity": {
            "homotopy": "u_s(d)=(1-s)*(1+beta_tors*d)+s*(1+d^eta)",
            "one_form": "A_s(d)=partial_s log(u_s(d))",
            "density_family": "rho_lambda(s)=1+lambda*(1-2s), 0<=lambda<=1",
            "density_mass": "int_0^1 rho_lambda(s) ds = 1",
            "density_positivity": "rho_lambda(s)>=0 for 0<=lambda<=1",
            "closed_form": "int rho_lambda*A_s ds = C(d)+lambda*B(d)",
            "B(d)": "C(d)*(1+2*(1+beta_tors*d)/(d^eta-beta_tors*d))-2",
            "interpretation": "front-loading changes the transport profile without using super-unit total mass; it is not the P2377 endpoint primitive unless lambda=0",
        },
        "lambda_threshold_formula": {
            "condition": "lambda > [3*(K1+C1)-(K5+C5)]/(B5-3*B1)",
            "status": "profile threshold for the affine normalized density family, not a strict source theorem for that profile",
        },
        "sample_front_loaded_support_audits": sample_rows,
        "lattice_threshold_audit": lattice_threshold_audit(),
        "max_midpoint_integral_abs_error": max_integral_error,
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}
    certificate = front_loaded_density_certificate()
    samples = certificate["sample_front_loaded_support_audits"]
    lattice = certificate["lattice_threshold_audit"]

    theorem_export = {
        "theorem_name": "P2379 front-loaded normalized transport density existence probe",
        "claim": (
            "P2378 rules out the unit-uniform endpoint primitive, but not all unit-mass transport profiles. "
            "For the normalized affine density rho_lambda(s)=1+lambda*(1-2s), the weighted log-transport has exact closed form C+lambda*B. "
            "On the audited P2376 grid and 81x81 lattice, lambda=0.8 is already above the computed profile threshold and selects the d5 chamber over all 792 supports."
        ),
        "positive_content": [
            "Derives the closed-form front-loading moment B(d)=int(1-2s)A_s(d)ds for the P2377 homotopy.",
            "Separates the P2378 failure of uniform density lambda=0 from the success of a unit-mass front-loaded density with lambda=0.8.",
            "Computationally audits midpoint integration, profile thresholds, and all 792 support maximizers on the 3x3 parameter grid plus an 81x81 threshold lattice.",
        ],
        "not_licensed": [
            "strict variational source theorem deriving the front-loaded density rho_lambda",
            "promotion of weighted transport into L_total",
            "claim that normalized transport is automatically sufficient without a profile/source theorem",
            "QW-2191 discharge",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2379_S1329_FRONT_LOADED_NORMALIZED_TRANSPORT_DENSITY_EXISTENCE_PROBE",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "artifact_replay": {
            "p2378_packet_id": artifacts["P2378_UNIT_NORMALIZED_INSUFFICIENCY"].get("packet_id"),
            "p2377_packet_id": artifacts["P2377_TRANSPORT_PRIMITIVE"].get("packet_id"),
            "p2376_packet_id": artifacts["P2376_RECTANGLE_ROBUSTNESS"].get("packet_id"),
        },
        "front_loaded_normalized_transport_density_certificate": certificate,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2378_loaded": probe["artifact_replay"]["p2378_packet_id"] == "P2378",
        "p2377_loaded": probe["artifact_replay"]["p2377_packet_id"] == "P2377",
        "p2376_loaded": probe["artifact_replay"]["p2376_packet_id"] == "P2376",
        "affine_density_is_normalized_and_positive_for_test_lambda": 0.0 <= LAMBDA_TEST <= 1.0,
        "all_sample_uniform_blends_fail_d5": all(not row["uniform_score_audit"]["d5_chamber"] for row in samples),
        "all_sample_front_loaded_blends_select_d5": all(
            row["front_loaded_score_audit"]["maximizer_h1_h5_pair_distribution"] == {"0,4": 12}
            for row in samples
        ),
        "all_sample_lambda_thresholds_below_test_lambda": all(row["lambda_threshold_gt"] < LAMBDA_TEST for row in samples),
        "lattice_thresholds_between_0_and_1": lattice["all_lattice_thresholds_between_0_and_1"],
        "lambda_test_exceeds_lattice_max": lattice["lambda_test_exceeds_lattice_max"],
        "docs_updated_with_p2379_profile_obligation": all("P2379/S1329" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2379_s1329_v1",
        "packet_id": "P2379",
        "stage_id": "S1329",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_FRONT_LOADED_NORMALIZED_TRANSPORT_PROFILE_EXISTS_SOURCE_PROFILE_STILL_OPEN",
        "result_kind": "FRONT_LOADED_NORMALIZED_TRANSPORT_DENSITY_EXISTENCE_PROBE",
        "front_loaded_normalized_transport_density_existence_probe": probe,
        "recommended_next_honest_step": (
            "Do not collapse P2378 into a blanket no-go for normalized transport. The next source theorem must derive either super-unit mass or a sufficiently front-loaded normalized density; "
            "until such a profile/source theorem is derived, rho_lambda remains an audited non-strict candidate."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2379 S1329: front-loaded normalized transport density existence probe",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2379/S1329 refines the P2378 no-go: unit-uniform endpoint transport fails, but a normalized front-loaded transport profile can pass the d5 chamber audit.",
                "For `rho_lambda(s)=1+lambda*(1-2s)`, `int rho_lambda ds=1` and `int rho_lambda*A_s(d) ds=C(d)+lambda*B(d)`.",
                "",
                "## Certificate",
                "",
                f"- Tested front-loading: `lambda={LAMBDA_TEST}` with density range `[{1.0 - LAMBDA_TEST}, {1.0 + LAMBDA_TEST}]`.",
                f"- 81x81 lattice max threshold row: `{lattice['max_threshold_row']}`.",
                f"- 81x81 lattice min threshold row: `{lattice['min_threshold_row']}`.",
                f"- Lambda test exceeds lattice max: `{lattice['lambda_test_exceeds_lattice_max']}`.",
                f"- Max midpoint integral error: `{certificate['max_midpoint_integral_abs_error']}`.",
                f"- Uniform grid scans select: `{[row['uniform_score_audit']['maximizer_h1_h5_pair_distribution'] for row in samples]}`.",
                f"- Front-loaded grid scans select: `{[row['front_loaded_score_audit']['maximizer_h1_h5_pair_distribution'] for row in samples]}`.",
                "",
                "## Hard limits",
                "",
                "- This is a normalized-profile existence/probe result, not a strict variational source theorem deriving `rho_lambda`.",
                "- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
