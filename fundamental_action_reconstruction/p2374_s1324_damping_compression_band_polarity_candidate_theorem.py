#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
from itertools import combinations
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2374_s1324_damping_compression_band_polarity_candidate_theorem.json"
MD = GEN / "p2374_s1324_damping_compression_band_polarity_candidate_theorem.md"

SOURCE_FILES = {
    "P2363_BRIDGE_MOMENT_TRANSPORT": GEN / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.json",
    "P2372_DIRECT_KERNEL_AUDIT": GEN / "p2372_s1322_bridge_kernel_direct_band_polarity_audit.json",
    "P2373_CORRECTION_CONE": GEN / "p2373_s1323_bridge_kernel_polarity_correction_cone_theorem.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
CHAMBER_THRESHOLD = 1.0 / 3.0
STRICT_PARAMS = {"omega": 743.0 / 4000.0, "phi": 13.0 / 80.0, "beta": 1.0, "eta": 9.0 / 5.0}
LEGACY_PARAMS = {"beta_tors": 1.0 / 100.0}
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


def compression_log_weight(d: int) -> float:
    """Strict nonlinear damping/compression surplus over legacy linear torsion damping."""
    strict_damping = 1.0 + STRICT_PARAMS["beta"] * d ** STRICT_PARAMS["eta"]
    legacy_damping = 1.0 + LEGACY_PARAMS["beta_tors"] * d
    return math.log(strict_damping / legacy_damping)


def internal_edges(support: tuple[int, ...], step: int) -> int:
    support_set = set(support)
    edges = set()
    for node in support_set:
        for neighbor in ((node + step) % Z12_NODE_COUNT, (node - step) % Z12_NODE_COUNT):
            if neighbor in support_set:
                edges.add(tuple(sorted((node, neighbor))))
    return len(edges)


def support_rows() -> list[dict[str, Any]]:
    rows = []
    for support in combinations(range(Z12_NODE_COUNT), SUPPORT_SIZE):
        rows.append({"support": list(support), "h1": internal_edges(support, 1), "h5": internal_edges(support, 5)})
    return rows


def score_maximizers(rows: list[dict[str, Any]], a: float, b: float) -> dict[str, Any]:
    scored = [(a * row["h1"] + b * row["h5"], row) for row in rows]
    maximum = max(score for score, _ in scored)
    maximizers = [row for score, row in scored if abs(score - maximum) <= 1e-10]
    pair_distribution: dict[str, int] = {}
    for row in maximizers:
        key = f"{row['h1']},{row['h5']}"
        pair_distribution[key] = pair_distribution.get(key, 0) + 1
    return {
        "weights": {"a_h1": a, "b_h5": b},
        "a_over_b": a / b if b > 0 else None,
        "d5_chamber": b > 0 and a >= 0 and a / b < CHAMBER_THRESHOLD,
        "maximum_score": maximum,
        "maximizer_count": len(maximizers),
        "maximizer_h1_h5_pair_distribution": dict(sorted(pair_distribution.items())),
        "maximizers_sample": maximizers[:12],
    }


def damping_compression_certificate() -> dict[str, Any]:
    rows = support_rows()
    k1 = k_strict(1)
    k5 = k_strict(5)
    c1 = compression_log_weight(1)
    c5 = compression_log_weight(5)
    denominator = c5 - 3.0 * c1
    tau_min = (3.0 * k1 - k5) / denominator
    tau = tau_min + EPSILON
    blended_a = k1 + tau * c1
    blended_b = k5 + tau * c5

    return {
        "strict_pair_weights": {"K1": k1, "K5": k5, "K1_over_K5": k1 / k5},
        "compression_log_weights": {
            "C1": c1,
            "C5": c5,
            "C1_over_C5": c1 / c5,
            "definition": "C(d)=log((1+beta*d^eta)/(1+beta_tors*d))",
        },
        "closed_form_chamber_tests": {
            "compression_alone_condition": "C1/C5 < 1/3",
            "compression_alone_condition_holds": c5 > 0 and c1 >= 0 and c1 / c5 < CHAMBER_THRESHOLD,
            "blend_condition": "(K1+tau*C1)/(K5+tau*C5)<1/3 iff tau>(3*K1-K5)/(C5-3*C1)",
            "blend_denominator_C5_minus_3C1": denominator,
            "blend_tau_gt": tau_min,
            "blend_tau_tested": tau,
        },
        "score_audits": {
            "strict_direct_kernel": score_maximizers(rows, k1, k5),
            "compression_log_weight_alone": score_maximizers(rows, c1, c5),
            "strict_plus_minimal_compression_blend_plus_epsilon": score_maximizers(rows, blended_a, blended_b),
        },
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}
    certificate = damping_compression_certificate()
    audits = certificate["score_audits"]

    theorem_export = {
        "theorem_name": "P2374 damping-compression band-polarity candidate theorem",
        "claim": (
            "The strict-side damping/compression surplus C(d)=log((1+beta*d^eta)/(1+beta_tors*d)) is, as a finite pair-feature, "
            "inside the P2371 d5 chamber at d=1,5. With the current canonical constants C1/C5 is about "
            f"{certificate['compression_log_weights']['C1_over_C5']:.6g}. Therefore C alone selects the 12 d5 path supports, and a positive blend "
            "K_strict(d)+tau*C(d) enters the d5 chamber exactly when "
            f"tau>{certificate['closed_form_chamber_tests']['blend_tau_gt']:.6g}. This is a candidate source direction, not a strict dynamical source theorem."
        ),
        "positive_content": [
            "Identifies the strict nonlinear damping/compression surplus as the first audited bridge-side feature with the correct d5 polarity.",
            "Derives the exact blend threshold against the P2372 direct strict-kernel failure.",
            "Verifies both compression-alone and just-above-threshold blended scores on all 792 five-node supports.",
        ],
        "not_licensed": [
            "strict dynamical source theorem for the compression feature",
            "promotion of C(d) into L_total or selector action without a source theorem",
            "QW-2191 discharge",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2374_S1324_DAMPING_COMPRESSION_BAND_POLARITY_CANDIDATE_THEOREM",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "artifact_replay": {
            "p2363_packet_id": artifacts["P2363_BRIDGE_MOMENT_TRANSPORT"].get("packet_id"),
            "p2372_packet_id": artifacts["P2372_DIRECT_KERNEL_AUDIT"].get("packet_id"),
            "p2373_packet_id": artifacts["P2373_CORRECTION_CONE"].get("packet_id"),
        },
        "damping_compression_certificate": certificate,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2363_loaded": probe["artifact_replay"]["p2363_packet_id"] == "P2363",
        "p2372_loaded": probe["artifact_replay"]["p2372_packet_id"] == "P2372",
        "p2373_loaded": probe["artifact_replay"]["p2373_packet_id"] == "P2373",
        "strict_direct_still_selects_h1_paths": audits["strict_direct_kernel"]["maximizer_h1_h5_pair_distribution"] == {"4,0": 12},
        "compression_log_weight_is_in_d5_chamber": certificate["closed_form_chamber_tests"]["compression_alone_condition_holds"],
        "compression_log_weight_selects_d5_paths": audits["compression_log_weight_alone"]["maximizer_h1_h5_pair_distribution"] == {"0,4": 12},
        "minimal_positive_blend_selects_d5_paths": audits["strict_plus_minimal_compression_blend_plus_epsilon"]["maximizer_h1_h5_pair_distribution"] == {"0,4": 12},
        "blend_threshold_is_finite_positive": certificate["closed_form_chamber_tests"]["blend_tau_gt"] > 0,
        "docs_updated_with_p2374_candidate": all("P2374/S1324" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2374_s1324_v1",
        "packet_id": "P2374",
        "stage_id": "S1324",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_DAMPING_COMPRESSION_POLARITY_CANDIDATE_NO_QW2191_DISCHARGE",
        "result_kind": "DAMPING_COMPRESSION_BAND_POLARITY_CANDIDATE_THEOREM",
        "damping_compression_band_polarity_candidate_theorem": probe,
        "recommended_next_honest_step": (
            "Try to export a variational or transport-level theorem that turns the audited compression surplus C(d) into an actual "
            "selector/source action. If that theorem is absent, keep C(d) as a quantified non-strict candidate rather than a closure."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2374 S1324: damping-compression band-polarity candidate theorem",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2374/S1324 audits the strict-side nonlinear damping/compression surplus as a candidate source direction for the P2373 polarity cone.",
                "",
                "## Certificate",
                "",
                f"- Compression definition: `{certificate['compression_log_weights']['definition']}`.",
                f"- Compression C1/C5: `{certificate['compression_log_weights']['C1_over_C5']}`; required `< 1/3`.",
                f"- Blend threshold tau must exceed: `{certificate['closed_form_chamber_tests']['blend_tau_gt']}`.",
                f"- Strict direct maximizers: `{audits['strict_direct_kernel']['maximizer_h1_h5_pair_distribution']}`.",
                f"- Compression-alone maximizers: `{audits['compression_log_weight_alone']['maximizer_h1_h5_pair_distribution']}`.",
                f"- Just-above-threshold blend maximizers: `{audits['strict_plus_minimal_compression_blend_plus_epsilon']['maximizer_h1_h5_pair_distribution']}`.",
                "",
                "## Hard limits",
                "",
                "- This is a bridge-side candidate direction and finite chamber theorem, not a strict dynamical source theorem.",
                "- No `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, L_total promotion, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
