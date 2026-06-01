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

OUT = GEN / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.json"
MD = GEN / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.md"

SOURCE_FILES = {
    "P2374_DAMPING_COMPRESSION_CANDIDATE": GEN / "p2374_s1324_damping_compression_band_polarity_candidate_theorem.json",
    "P2375_INTERVAL_ROBUSTNESS": GEN / "p2375_s1325_damping_compression_polarity_interval_robustness_theorem.json",
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
EPSILON = 1e-9
QUADRATURE_STEPS = 80000


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


def transport_speed(d: int, eta: float, beta_tors: float, s: float) -> float:
    """Exact log-speed for the linear denominator homotopy from legacy to strict damping."""
    legacy = 1.0 + beta_tors * d
    strict = 1.0 + STRICT_BETA * d**eta
    return (strict - legacy) / ((1.0 - s) * legacy + s * strict)


def midpoint_integral_transport_speed(d: int, eta: float, beta_tors: float, steps: int = QUADRATURE_STEPS) -> float:
    total = 0.0
    for index in range(steps):
        s = (index + 0.5) / steps
        total += transport_speed(d, eta, beta_tors, s)
    return total / steps


def blend_denominator(eta: float, beta_tors: float) -> float:
    return compression_log_weight(5, eta, beta_tors) - 3.0 * compression_log_weight(1, eta, beta_tors)


def denominator_dx(eta: float, beta_tors: float) -> float:
    # d/dx [log((1+5^eta)/(1+5x)) - 3 log(2/(1+x))]
    return -5.0 / (1.0 + 5.0 * beta_tors) + 3.0 / (1.0 + beta_tors)


def denominator_deta(eta: float, beta_tors: float) -> float:
    return math.log(5.0) * 5.0**eta / (1.0 + 5.0**eta)


def uniform_tau_threshold() -> dict[str, Any]:
    eta_min, eta_max = ETA_INTERVAL
    beta_min, beta_max = BETA_TORS_INTERVAL
    k1 = k_strict(1)
    k5 = k_strict(5)
    numerator = 3.0 * k1 - k5
    denominator_corner = blend_denominator(eta_min, beta_max)
    return {
        "strict_pair_weights": {"K1": k1, "K5": k5, "3K1_minus_K5": numerator, "K1_over_K5": k1 / k5},
        "denominator": {
            "D(eta,x)": "C5(eta,x)-3*C1(eta,x)",
            "dD_deta": "ln(5)*5^eta/(1+5^eta) > 0",
            "dD_dx": "-5/(1+5*x)+3/(1+x) <= 0 on x in [0,0.1]",
            "minimum_corner": {"eta": eta_min, "beta_tors": beta_max},
            "minimum_corner_value": denominator_corner,
            "minimum_corner_dD_deta": denominator_deta(eta_min, beta_max),
            "right_endpoint_dD_dx": denominator_dx(eta_min, beta_max),
            "left_endpoint_dD_dx": denominator_dx(eta_min, beta_min),
            "positive_on_rectangle": denominator_corner > 0,
        },
        "uniform_coupling": {
            "condition": "tau > (3*K1-K5)/min_rectangle(C5-3*C1)",
            "tau_gt_uniform": numerator / denominator_corner,
            "tau_tested": numerator / denominator_corner + EPSILON,
            "normalization_status": "scalar coupling still not derived by strict dynamics",
        },
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


def transport_certificate() -> dict[str, Any]:
    rows = support_rows()
    threshold = uniform_tau_threshold()
    tau = threshold["uniform_coupling"]["tau_tested"]
    k1 = threshold["strict_pair_weights"]["K1"]
    k5 = threshold["strict_pair_weights"]["K5"]

    samples = []
    max_integral_error = 0.0
    for eta, beta_tors in product(SAMPLE_ETAS, SAMPLE_BETA_TORS):
        c1 = compression_log_weight(1, eta, beta_tors)
        c5 = compression_log_weight(5, eta, beta_tors)
        i1 = midpoint_integral_transport_speed(1, eta, beta_tors)
        i5 = midpoint_integral_transport_speed(5, eta, beta_tors)
        max_integral_error = max(max_integral_error, abs(i1 - c1), abs(i5 - c5))
        blended_a = k1 + tau * c1
        blended_b = k5 + tau * c5
        samples.append(
            {
                "eta": eta,
                "beta_tors": beta_tors,
                "C1_endpoint_log": c1,
                "C5_endpoint_log": c5,
                "C1_midpoint_integral": i1,
                "C5_midpoint_integral": i5,
                "C1_integral_abs_error": abs(i1 - c1),
                "C5_integral_abs_error": abs(i5 - c5),
                "blend_denominator_C5_minus_3C1": c5 - 3.0 * c1,
                "blended_score_audit": score_maximizers(rows, blended_a, blended_b),
            }
        )

    return {
        "transport_identity": {
            "legacy_damping_endpoint": "u_0(d)=1+beta_tors*d",
            "strict_damping_endpoint": "u_1(d)=1+d^eta",
            "homotopy": "u_s(d)=(1-s)*(1+beta_tors*d)+s*(1+d^eta)",
            "one_form": "A_s(d)=partial_s log(u_s(d))=((1+d^eta)-(1+beta_tors*d))/u_s(d)",
            "primitive": "integral_0^1 A_s(d) ds = log((1+d^eta)/(1+beta_tors*d)) = C(d)",
            "path_status": "exact endpoint transport primitive for damping completion, not a variational source by itself",
        },
        "uniform_threshold_certificate": threshold,
        "sample_transport_and_support_audits": samples,
        "max_midpoint_integral_abs_error": max_integral_error,
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}
    certificate = transport_certificate()
    samples = certificate["sample_transport_and_support_audits"]
    threshold = certificate["uniform_threshold_certificate"]

    theorem_export = {
        "theorem_name": "P2377 damping-compression transport primitive and uniform coupling theorem",
        "claim": (
            "The audited compression feature C(d) is the exact endpoint primitive of the damping-completion log-transport one-form "
            "A_s(d)=partial_s log((1-s)*(1+beta_tors*d)+s*(1+d^eta)). On the P2376 rectangle, the blend denominator "
            "C5-3*C1 is minimized at eta=9/5, beta_tors=0.1, giving a uniform scalar-coupling threshold for overcoming the failed direct strict-kernel pair weights."
        ),
        "positive_content": [
            "Exports a closed transport-primitive provenance for C(d), instead of treating C(d) as a numerological pair feature.",
            "Computes a rectangle-uniform coupling acceptance threshold for K(d)+tau*C(d).",
            "Checks endpoint-integral equality and d5 support selection on a 3x3 eta/beta_tors grid over all 792 supports.",
        ],
        "not_licensed": [
            "strict variational source theorem fixing the scalar coupling tau",
            "promotion of C(d) or tau*C(d) into L_total without a separate variational/source theorem",
            "derivation of eta, beta, beta_tors, omega, or phi from strict nadsoliton dynamics",
            "QW-2191 discharge",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2377_S1327_DAMPING_COMPRESSION_TRANSPORT_PRIMITIVE_UNIFORM_COUPLING_THEOREM",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "artifact_replay": {
            "p2374_packet_id": artifacts["P2374_DAMPING_COMPRESSION_CANDIDATE"].get("packet_id"),
            "p2375_packet_id": artifacts["P2375_INTERVAL_ROBUSTNESS"].get("packet_id"),
            "p2376_packet_id": artifacts["P2376_RECTANGLE_ROBUSTNESS"].get("packet_id"),
        },
        "transport_primitive_uniform_coupling_certificate": certificate,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2374_loaded": probe["artifact_replay"]["p2374_packet_id"] == "P2374",
        "p2375_loaded": probe["artifact_replay"]["p2375_packet_id"] == "P2375",
        "p2376_loaded": probe["artifact_replay"]["p2376_packet_id"] == "P2376",
        "transport_integral_matches_endpoint_log_on_grid": certificate["max_midpoint_integral_abs_error"] < 1e-8,
        "uniform_denominator_positive": threshold["denominator"]["positive_on_rectangle"],
        "uniform_tau_positive": threshold["uniform_coupling"]["tau_gt_uniform"] > 0,
        "all_uniform_blend_scans_select_d5_paths": all(
            row["blended_score_audit"]["maximizer_h1_h5_pair_distribution"] == {"0,4": 12} for row in samples
        ),
        "docs_updated_with_p2377_transport": all("P2377/S1327" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2377_s1327_v1",
        "packet_id": "P2377",
        "stage_id": "S1327",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_TRANSPORT_PRIMITIVE_AND_UNIFORM_COUPLING_NO_VARIATIONAL_SOURCE_CLOSURE",
        "result_kind": "DAMPING_COMPRESSION_TRANSPORT_PRIMITIVE_UNIFORM_COUPLING_THEOREM",
        "damping_compression_transport_primitive_uniform_coupling_theorem": probe,
        "recommended_next_honest_step": (
            "Use P2377 only as transport provenance plus a uniform acceptance threshold: next either derive a genuine variational/source theorem "
            "that fixes a coupling tau above the threshold, or mark tau*C(d) as an explicit non-strict selector premise."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2377 S1327: damping-compression transport primitive and uniform coupling theorem",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2377/S1327 proves that `C(d)=log((1+d^eta)/(1+beta_tors*d))` is the exact endpoint primitive of the damping-completion log-transport one-form.",
                "It also computes a rectangle-uniform coupling threshold for the blended pair score `K_strict(d)+tau*C(d)`.",
                "",
                "## Certificate",
                "",
                f"- Transport primitive: `{certificate['transport_identity']['primitive']}`.",
                f"- Uniform denominator corner: `{threshold['denominator']['minimum_corner']}`.",
                f"- Uniform denominator value: `{threshold['denominator']['minimum_corner_value']}`.",
                f"- Uniform tau threshold: `{threshold['uniform_coupling']['tau_gt_uniform']}`.",
                f"- Max midpoint integral error on grid: `{certificate['max_midpoint_integral_abs_error']}`.",
                f"- Grid blend scans select: `{[row['blended_score_audit']['maximizer_h1_h5_pair_distribution'] for row in samples]}`.",
                "",
                "## Hard limits",
                "",
                "- This is transport provenance plus a uniform coupling acceptance theorem, not a strict variational source theorem fixing `tau`.",
                "- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
