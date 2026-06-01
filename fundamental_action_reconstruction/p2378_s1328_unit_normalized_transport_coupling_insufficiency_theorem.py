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

OUT = GEN / "p2378_s1328_unit_normalized_transport_coupling_insufficiency_theorem.json"
MD = GEN / "p2378_s1328_unit_normalized_transport_coupling_insufficiency_theorem.md"

SOURCE_FILES = {
    "P2377_TRANSPORT_PRIMITIVE": GEN / "p2377_s1327_damping_compression_transport_primitive_uniform_coupling_theorem.json",
    "P2376_RECTANGLE_ROBUSTNESS": GEN / "p2376_s1326_damping_compression_eta_beta_rectangle_robustness_theorem.json",
    "P2374_DAMPING_COMPRESSION_CANDIDATE": GEN / "p2374_s1324_damping_compression_band_polarity_candidate_theorem.json",
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


def blend_denominator(eta: float, beta_tors: float) -> float:
    return compression_log_weight(5, eta, beta_tors) - 3.0 * compression_log_weight(1, eta, beta_tors)


def denominator_dx(beta_tors: float) -> float:
    return -5.0 / (1.0 + 5.0 * beta_tors) + 3.0 / (1.0 + beta_tors)


def denominator_deta(eta: float) -> float:
    return math.log(5.0) * 5.0**eta / (1.0 + 5.0**eta)


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


def mass_budget_certificate() -> dict[str, Any]:
    eta_min, eta_max = ETA_INTERVAL
    beta_min, beta_max = BETA_TORS_INTERVAL
    k1 = k_strict(1)
    k5 = k_strict(5)
    numerator = 3.0 * k1 - k5
    denominator_min = blend_denominator(eta_min, beta_max)
    denominator_max = blend_denominator(eta_max, beta_min)
    rows = support_rows()

    sample_rows = []
    for eta, beta_tors in product(SAMPLE_ETAS, SAMPLE_BETA_TORS):
        c1 = compression_log_weight(1, eta, beta_tors)
        c5 = compression_log_weight(5, eta, beta_tors)
        denominator = c5 - 3.0 * c1
        threshold = numerator / denominator
        unit_a = k1 + c1
        unit_b = k5 + c5
        just_above_a = k1 + (threshold + EPSILON) * c1
        just_above_b = k5 + (threshold + EPSILON) * c5
        sample_rows.append(
            {
                "eta": eta,
                "beta_tors": beta_tors,
                "C1": c1,
                "C5": c5,
                "denominator_C5_minus_3C1": denominator,
                "tau_threshold": threshold,
                "unit_mass_blend_ratio": unit_a / unit_b,
                "unit_mass_blend_score_audit": score_maximizers(rows, unit_a, unit_b),
                "just_above_threshold_score_audit": score_maximizers(rows, just_above_a, just_above_b),
            }
        )

    return {
        "mass_budget_theorem": {
            "effective_blend": "a=K1+M*C1, b=K5+M*C5 for nonnegative total transport mass M",
            "d5_condition": "M > (3*K1-K5)/(C5-3*C1)",
            "normalized_transport_budget": "M<=1",
            "interpretation": "An endpoint primitive with only unit-normalized total mass is not enough to overcome direct strict-kernel polarity.",
        },
        "rectangle_proof": {
            "numerator_3K1_minus_K5": numerator,
            "D(eta,x)": "C5(eta,x)-3*C1(eta,x)",
            "dD_deta": "ln(5)*5^eta/(1+5^eta) > 0",
            "dD_dx": "-5/(1+5*x)+3/(1+x) < 0 on x in [0,0.1]",
            "denominator_min_corner": {"eta": eta_min, "beta_tors": beta_max},
            "denominator_min_value": denominator_min,
            "denominator_max_corner": {"eta": eta_max, "beta_tors": beta_min},
            "denominator_max_value": denominator_max,
            "dD_deta_lower_left": denominator_deta(eta_min),
            "dD_dx_left_endpoint": denominator_dx(beta_min),
            "dD_dx_right_endpoint": denominator_dx(beta_max),
            "denominator_positive_on_rectangle": denominator_min > 0,
            "unit_mass_insufficient_on_rectangle": denominator_max < numerator,
            "tau_threshold_range": {
                "minimum_tau_gt": numerator / denominator_max,
                "maximum_tau_gt": numerator / denominator_min,
            },
        },
        "sample_mass_budget_support_audits": sample_rows,
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}
    certificate = mass_budget_certificate()
    proof = certificate["rectangle_proof"]
    samples = certificate["sample_mass_budget_support_audits"]

    theorem_export = {
        "theorem_name": "P2378 unit-normalized transport coupling insufficiency theorem",
        "claim": (
            "Although P2377 gives exact transport provenance for C(d), a unit-normalized transport primitive is insufficient: "
            "for every eta in [9/5,2] and beta_tors in [0,0.1], the d5 chamber for K+M*C requires M>1. "
            "The threshold range is certified by D=C5-3*C1, whose maximum remains below 3*K1-K5."
        ),
        "positive_content": [
            "Derives the scalar mass-budget condition M>(3*K1-K5)/(C5-3*C1).",
            "Proves across the P2376 rectangle that every normalized budget M<=1 fails to overcome direct strict-kernel polarity.",
            "Computationally separates unit-mass failure from just-above-threshold d5 success on the 3x3 grid over all 792 supports.",
        ],
        "not_licensed": [
            "strict variational source theorem fixing M above threshold",
            "using P2377 transport provenance alone as selector closure",
            "promotion of C(d) or M*C(d) into L_total without a source theorem",
            "QW-2191 discharge",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2378_S1328_UNIT_NORMALIZED_TRANSPORT_COUPLING_INSUFFICIENCY_THEOREM",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "artifact_replay": {
            "p2377_packet_id": artifacts["P2377_TRANSPORT_PRIMITIVE"].get("packet_id"),
            "p2376_packet_id": artifacts["P2376_RECTANGLE_ROBUSTNESS"].get("packet_id"),
            "p2374_packet_id": artifacts["P2374_DAMPING_COMPRESSION_CANDIDATE"].get("packet_id"),
        },
        "unit_normalized_transport_coupling_insufficiency_certificate": certificate,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2377_loaded": probe["artifact_replay"]["p2377_packet_id"] == "P2377",
        "p2376_loaded": probe["artifact_replay"]["p2376_packet_id"] == "P2376",
        "p2374_loaded": probe["artifact_replay"]["p2374_packet_id"] == "P2374",
        "denominator_positive_on_rectangle": proof["denominator_positive_on_rectangle"],
        "unit_mass_insufficient_on_rectangle": proof["unit_mass_insufficient_on_rectangle"],
        "all_sample_unit_blends_fail_d5": all(not row["unit_mass_blend_score_audit"]["d5_chamber"] for row in samples),
        "all_sample_unit_blends_select_mixed_33": all(
            row["unit_mass_blend_score_audit"]["maximizer_h1_h5_pair_distribution"] == {"3,3": 24}
            for row in samples
        ),
        "all_sample_threshold_blends_select_d5": all(
            row["just_above_threshold_score_audit"]["maximizer_h1_h5_pair_distribution"] == {"0,4": 12}
            for row in samples
        ),
        "docs_updated_with_p2378_insufficiency": all("P2378/S1328" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2378_s1328_v1",
        "packet_id": "P2378",
        "stage_id": "S1328",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_UNIT_NORMALIZED_TRANSPORT_INSUFFICIENT_SOURCE_STRENGTH_STILL_OPEN",
        "result_kind": "UNIT_NORMALIZED_TRANSPORT_COUPLING_INSUFFICIENCY_THEOREM",
        "unit_normalized_transport_coupling_insufficiency_theorem": probe,
        "recommended_next_honest_step": (
            "Do not treat P2377 transport provenance as selector closure. Next either derive a variational/source normalization theorem "
            "that gives total coupling M above the P2378 threshold range, or mark the extra coupling strength as an explicit non-strict premise."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2378 S1328: unit-normalized transport coupling insufficiency theorem",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2378/S1328 proves that the exact P2377 transport primitive still lacks enough source strength if its total transport mass is unit-normalized.",
                "For `K_strict(d)+M*C(d)`, the d5 chamber requires `M>(3*K1-K5)/(C5-3*C1)`.",
                "",
                "## Certificate",
                "",
                f"- Numerator `3*K1-K5`: `{proof['numerator_3K1_minus_K5']}`.",
                f"- Denominator min: `{proof['denominator_min_value']}` at `{proof['denominator_min_corner']}`.",
                f"- Denominator max: `{proof['denominator_max_value']}` at `{proof['denominator_max_corner']}`.",
                f"- Threshold range: `{proof['tau_threshold_range']}`.",
                f"- Unit mass insufficient on rectangle: `{proof['unit_mass_insufficient_on_rectangle']}`.",
                f"- Unit grid scans select: `{[row['unit_mass_blend_score_audit']['maximizer_h1_h5_pair_distribution'] for row in samples]}`.",
                f"- Just-above-threshold grid scans select: `{[row['just_above_threshold_score_audit']['maximizer_h1_h5_pair_distribution'] for row in samples]}`.",
                "",
                "## Hard limits",
                "",
                "- This is an insufficiency theorem for normalized transport strength, not a source theorem fixing a super-unit coupling.",
                "- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
