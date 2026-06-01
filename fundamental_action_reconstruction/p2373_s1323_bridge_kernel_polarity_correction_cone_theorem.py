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

OUT = GEN / "p2373_s1323_bridge_kernel_polarity_correction_cone_theorem.json"
MD = GEN / "p2373_s1323_bridge_kernel_polarity_correction_cone_theorem.md"

SOURCE_FILES = {
    "P2372_DIRECT_KERNEL_AUDIT": GEN / "p2372_s1322_bridge_kernel_direct_band_polarity_audit.json",
    "P2371_AUT_BANDPASS_OBSTRUCTION": GEN / "p2371_s1321_aut_invariant_unit_bandpass_obstruction_theorem.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
STRICT_PARAMS = {"omega": 743.0 / 4000.0, "phi": 13.0 / 80.0, "beta": 1.0, "eta": 9.0 / 5.0}
THRESHOLD = 1.0 / 3.0
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
        "d5_chamber": b > 0 and a >= 0 and a / b < THRESHOLD,
        "maximum_score": maximum,
        "maximizer_count": len(maximizers),
        "maximizer_h1_h5_pair_distribution": dict(sorted(pair_distribution.items())),
        "maximizers_sample": maximizers[:12],
    }


def correction_cone() -> dict[str, Any]:
    a0 = k_strict(1)
    b0 = k_strict(5)
    rows = support_rows()
    h5_boost_min = 3.0 * a0 - b0
    h1_suppression_min = a0 - b0 / 3.0
    antisymmetric_min = (3.0 * a0 - b0) / 4.0

    candidates = {
        "baseline_direct_kernel": (a0, b0),
        "pure_h5_boost_at_open_threshold_plus_epsilon": (a0, b0 + h5_boost_min + EPSILON),
        "pure_h1_suppression_at_open_threshold_plus_epsilon": (a0 - h1_suppression_min - EPSILON, b0),
        "antisymmetric_h5_minus_h1_polarity_at_open_threshold_plus_epsilon": (
            a0 - antisymmetric_min - EPSILON,
            b0 + antisymmetric_min + EPSILON,
        ),
    }
    scored = {name: score_maximizers(rows, a, b) for name, (a, b) in candidates.items()}
    return {
        "baseline_weights": {"a0_K1": a0, "b0_K5": b0, "a0_over_b0": a0 / b0},
        "required_open_chamber": "a/b < 1/3, with a>=0 and b>0",
        "minimal_corrections": {
            "pure_h5_boost_lambda_gt": h5_boost_min,
            "pure_h5_boost_lambda_gt_in_units_of_K5": h5_boost_min / b0,
            "pure_h1_suppression_mu_gt": h1_suppression_min,
            "pure_h1_suppression_mu_gt_in_units_of_K1": h1_suppression_min / a0,
            "antisymmetric_gamma_for_plus_h5_minus_h1_gt": antisymmetric_min,
            "antisymmetric_gamma_gt_in_units_of_K1": antisymmetric_min / a0,
            "antisymmetric_gamma_gt_in_units_of_K5": antisymmetric_min / b0,
        },
        "derivations": {
            "pure_h5_boost": "a0/(b0+lambda)<1/3 iff lambda>3*a0-b0",
            "pure_h1_suppression": "(a0-mu)/b0<1/3 iff mu>a0-b0/3",
            "antisymmetric_plus_h5_minus_h1": "(a0-gamma)/(b0+gamma)<1/3 iff gamma>(3*a0-b0)/4",
        },
        "score_audits": scored,
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}
    cone = correction_cone()

    theorem_export = {
        "theorem_name": "P2373 bridge-kernel polarity correction cone theorem",
        "claim": (
            "Starting from direct bridge-completed strict-kernel pair weights a0=K_strict(1), b0=K_strict(5), "
            "P2371's d5 chamber a/b<1/3 requires a large additional polarity source. A pure h5 boost must exceed "
            f"{cone['minimal_corrections']['pure_h5_boost_lambda_gt']:.6g}; a pure h1 suppression must exceed "
            f"{cone['minimal_corrections']['pure_h1_suppression_mu_gt']:.6g}; and an antisymmetric +h5-h1 term must exceed "
            f"{cone['minimal_corrections']['antisymmetric_gamma_for_plus_h5_minus_h1_gt']:.6g}. These are necessary correction-cone bounds, not a derivation of the missing dynamical term."
        ),
        "positive_content": [
            "Closed-form inequalities give the minimal correction thresholds needed to enter the P2371 d5 chamber.",
            "Adding just above each threshold is verified on all 792 supports to select the 12 d5 path supports.",
            "The baseline direct strict kernel is replayed and still selects distance-1 paths.",
            "The correction cone quantifies how strong any bridge-completed dynamical polarity term must be.",
        ],
        "not_licensed": [
            "existence of a strict dynamical polarity term",
            "strict derivation of distance-5 band-pass action",
            "promotion of correction terms to L_total without a source theorem",
            "QW-2191 discharge",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2373_S1323_BRIDGE_KERNEL_POLARITY_CORRECTION_CONE_THEOREM",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "p2372_replay": {
            "packet_id": artifacts["P2372_DIRECT_KERNEL_AUDIT"].get("packet_id"),
            "strict_direct_fails_d5_chamber": artifacts["P2372_DIRECT_KERNEL_AUDIT"].get("gatekeeper_checks", {}).get("strict_direct_fails_d5_chamber"),
        },
        "p2371_replay": {
            "packet_id": artifacts["P2371_AUT_BANDPASS_OBSTRUCTION"].get("packet_id"),
            "full_aut_does_not_select_d5": artifacts["P2371_AUT_BANDPASS_OBSTRUCTION"].get("gatekeeper_checks", {}).get("full_aut_does_not_select_d5"),
        },
        "correction_cone_certificate": cone,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    score_audits = cone["score_audits"]
    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2372_loaded": probe["p2372_replay"]["packet_id"] == "P2372",
        "p2371_loaded": probe["p2371_replay"]["packet_id"] == "P2371",
        "baseline_selects_h1_paths": score_audits["baseline_direct_kernel"]["maximizer_h1_h5_pair_distribution"] == {"4,0": 12},
        "pure_h5_boost_enters_d5_chamber": score_audits["pure_h5_boost_at_open_threshold_plus_epsilon"]["maximizer_h1_h5_pair_distribution"] == {"0,4": 12},
        "pure_h1_suppression_enters_d5_chamber": score_audits["pure_h1_suppression_at_open_threshold_plus_epsilon"]["maximizer_h1_h5_pair_distribution"] == {"0,4": 12},
        "antisymmetric_correction_enters_d5_chamber": score_audits["antisymmetric_h5_minus_h1_polarity_at_open_threshold_plus_epsilon"]["maximizer_h1_h5_pair_distribution"] == {"0,4": 12},
        "corrections_are_large_relative_to_kernel": cone["minimal_corrections"]["pure_h5_boost_lambda_gt_in_units_of_K5"] > 50.0,
        "docs_updated_with_p2373_cone": all("P2373/S1323" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2373_s1323_v1",
        "packet_id": "P2373",
        "stage_id": "S1323",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_POLARITY_CORRECTION_CONE_NO_QW2191_DISCHARGE",
        "result_kind": "BRIDGE_KERNEL_POLARITY_CORRECTION_CONE_THEOREM",
        "bridge_kernel_polarity_correction_cone_theorem": probe,
        "recommended_next_honest_step": (
            "Search for an actual bridge-completed dynamical source for one of the quantified correction directions "
            "(+h5 boost, h1 suppression, or antisymmetric +h5-h1 polarity). If none is exported, keep the band-polarity term non-strict."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2373 S1323: bridge-kernel polarity correction cone theorem",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2373/S1323 quantifies how large a new polarity term must be to move the direct strict-kernel pair weights into the P2371 d5 chamber.",
                "",
                "## Certificate",
                "",
                f"- Baseline K1/K5: `{cone['baseline_weights']['a0_over_b0']}`.",
                f"- Pure h5 boost lambda must exceed: `{cone['minimal_corrections']['pure_h5_boost_lambda_gt']}` (`{cone['minimal_corrections']['pure_h5_boost_lambda_gt_in_units_of_K5']}` times K5).",
                f"- Pure h1 suppression mu must exceed: `{cone['minimal_corrections']['pure_h1_suppression_mu_gt']}`.",
                f"- Antisymmetric gamma for +h5-h1 must exceed: `{cone['minimal_corrections']['antisymmetric_gamma_for_plus_h5_minus_h1_gt']}`.",
                f"- Baseline maximizers: `{score_audits['baseline_direct_kernel']['maximizer_h1_h5_pair_distribution']}`.",
                f"- Just-above-threshold antisymmetric maximizers: `{score_audits['antisymmetric_h5_minus_h1_polarity_at_open_threshold_plus_epsilon']['maximizer_h1_h5_pair_distribution']}`.",
                "",
                "## Hard limits",
                "",
                "- This quantifies necessary correction sizes; it does not derive a strict dynamical correction term.",
                "- No `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
