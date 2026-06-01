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

OUT = GEN / "p2376_s1326_damping_compression_eta_beta_rectangle_robustness_theorem.json"
MD = GEN / "p2376_s1326_damping_compression_eta_beta_rectangle_robustness_theorem.md"

SOURCE_FILES = {
    "P2375_INTERVAL_ROBUSTNESS": GEN / "p2375_s1325_damping_compression_polarity_interval_robustness_theorem.json",
    "P2374_DAMPING_COMPRESSION_CANDIDATE": GEN / "p2374_s1324_damping_compression_band_polarity_candidate_theorem.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
ETA_INTERVAL = (9.0 / 5.0, 2.0)
BETA_TORS_INTERVAL = (0.0, 0.1)
STRICT_BETA = 1.0
CHAMBER_THRESHOLD = 1.0 / 3.0
SAMPLE_ETAS = [9.0 / 5.0, 19.0 / 10.0, 2.0]
SAMPLE_BETA_TORS = [0.0, 0.01, 0.1]


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


def compression_log_weight(d: int, eta: float, beta_tors: float) -> float:
    return math.log((1.0 + STRICT_BETA * d**eta) / (1.0 + beta_tors * d))


def margin(eta: float, beta_tors: float) -> float:
    return (1.0 + 5.0**eta) * (1.0 + beta_tors) ** 3 - 8.0 * (1.0 + 5.0 * beta_tors)


def margin_dx(eta: float, beta_tors: float) -> float:
    return 3.0 * (1.0 + 5.0**eta) * (1.0 + beta_tors) ** 2 - 40.0


def margin_deta(eta: float, beta_tors: float) -> float:
    return math.log(5.0) * 5.0**eta * (1.0 + beta_tors) ** 3


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
        "maximizer_count": len(maximizers),
        "maximizer_h1_h5_pair_distribution": dict(sorted(pair_distribution.items())),
    }


def rectangle_certificate() -> dict[str, Any]:
    eta_min, eta_max = ETA_INTERVAL
    beta_min, beta_max = BETA_TORS_INTERVAL
    rows = support_rows()
    sample_rows = []
    for eta, beta_tors in product(SAMPLE_ETAS, SAMPLE_BETA_TORS):
        c1 = compression_log_weight(1, eta, beta_tors)
        c5 = compression_log_weight(5, eta, beta_tors)
        sample_rows.append(
            {
                "eta": eta,
                "beta_tors": beta_tors,
                "C1": c1,
                "C5": c5,
                "C1_over_C5": c1 / c5,
                "margin": margin(eta, beta_tors),
                "score_audit": score_maximizers(rows, c1, c5),
            }
        )

    corner_margin = margin(eta_min, beta_min)
    corner_dx = margin_dx(eta_min, beta_min)
    corner_deta = margin_deta(eta_min, beta_min)
    return {
        "rectangle": {
            "eta_min": eta_min,
            "eta_max": eta_max,
            "beta_tors_min": beta_min,
            "beta_tors_max": beta_max,
        },
        "symbolic_inequality": {
            "target": "C(1;eta,x)/C(5;eta,x) < 1/3",
            "equivalent_positive_margin": "F(eta,x)=(1+5^eta)*(1+x)^3 - 8*(1+5*x) > 0",
            "dF_dx": "3*(1+5^eta)*(1+x)^2 - 40",
            "dF_deta": "ln(5)*5^eta*(1+x)^3",
            "minimum_corner": "eta=9/5, x=0",
            "minimum_corner_margin": corner_margin,
            "minimum_corner_dF_dx": corner_dx,
            "minimum_corner_dF_deta": corner_deta,
        },
        "proof_summary": {
            "dF_deta_positive_on_rectangle": corner_deta > 0,
            "dF_dx_positive_on_rectangle": corner_dx > 0,
            "margin_minimum_at_lower_left_corner": corner_deta > 0 and corner_dx > 0,
            "margin_positive_on_rectangle": corner_margin > 0 and corner_deta > 0 and corner_dx > 0,
            "therefore_d5_chamber_on_rectangle": corner_margin > 0 and corner_deta > 0 and corner_dx > 0,
        },
        "sample_score_audits": sample_rows,
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}
    certificate = rectangle_certificate()
    sample_rows = certificate["sample_score_audits"]

    theorem_export = {
        "theorem_name": "P2376 damping-compression eta-beta rectangle robustness theorem",
        "claim": (
            "For every eta in [9/5,2] and beta_tors in [0,0.1], the compression surplus "
            "C(d)=log((1+d^eta)/(1+beta_tors*d)) remains in the d5 chamber C(1)/C(5)<1/3. "
            "The positive-margin function F(eta,x)=(1+5^eta)*(1+x)^3-8*(1+5*x) is increasing in both eta and x on the rectangle, "
            "and its lower-left value is positive. This makes the P2374/P2375 candidate robust to the audited eta and beta_tors rectangle, without deriving a source theorem."
        ),
        "positive_content": [
            "Closed-form bivariate rectangle proof for the compression-polarity chamber inequality.",
            "Support scans over eta/beta_tors grid points verify d5 selection on all 792 supports.",
            "Strengthens the candidate against eta and beta_tors fine-tuning concerns.",
        ],
        "not_licensed": [
            "strict dynamical source theorem for C(d)",
            "derivation of eta=9/5, beta=1, or beta_tors from nadsoliton dynamics",
            "promotion of C(d) into L_total or selector action without a variational/transport theorem",
            "QW-2191 discharge",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2376_S1326_DAMPING_COMPRESSION_ETA_BETA_RECTANGLE_ROBUSTNESS_THEOREM",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "artifact_replay": {
            "p2375_packet_id": artifacts["P2375_INTERVAL_ROBUSTNESS"].get("packet_id"),
            "p2374_packet_id": artifacts["P2374_DAMPING_COMPRESSION_CANDIDATE"].get("packet_id"),
        },
        "rectangle_robustness_certificate": certificate,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2375_loaded": probe["artifact_replay"]["p2375_packet_id"] == "P2375",
        "p2374_loaded": probe["artifact_replay"]["p2374_packet_id"] == "P2374",
        "margin_positive_on_rectangle": certificate["proof_summary"]["margin_positive_on_rectangle"],
        "all_sample_ratios_in_d5_chamber": all(row["C1_over_C5"] < CHAMBER_THRESHOLD for row in sample_rows),
        "all_sample_scans_select_d5_paths": all(
            row["score_audit"]["maximizer_h1_h5_pair_distribution"] == {"0,4": 12} for row in sample_rows
        ),
        "docs_updated_with_p2376_rectangle": all("P2376/S1326" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2376_s1326_v1",
        "packet_id": "P2376",
        "stage_id": "S1326",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_DAMPING_COMPRESSION_ETA_BETA_RECTANGLE_ROBUST_NO_QW2191_DISCHARGE",
        "result_kind": "DAMPING_COMPRESSION_ETA_BETA_RECTANGLE_ROBUSTNESS_THEOREM",
        "damping_compression_eta_beta_rectangle_robustness_theorem": probe,
        "recommended_next_honest_step": (
            "Stop widening finite robustness unless a new parameter risk appears; next try to derive a variational/transport source theorem "
            "that licenses the robust C(d) polarity feature as an actual selector action, or keep it non-strict."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2376 S1326: damping-compression eta-beta rectangle robustness theorem",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2376/S1326 proves that the P2374/P2375 compression-polarity candidate is robust on `eta in [9/5,2]` and `beta_tors in [0,0.1]`.",
                "",
                "## Certificate",
                "",
                f"- Margin: `{certificate['symbolic_inequality']['equivalent_positive_margin']}`.",
                f"- Minimum corner margin: `{certificate['symbolic_inequality']['minimum_corner_margin']}`.",
                f"- Minimum corner dF/dx: `{certificate['symbolic_inequality']['minimum_corner_dF_dx']}`.",
                f"- Minimum corner dF/deta: `{certificate['symbolic_inequality']['minimum_corner_dF_deta']}`.",
                f"- Grid support scans select: `{[row['score_audit']['maximizer_h1_h5_pair_distribution'] for row in sample_rows]}`.",
                "",
                "## Hard limits",
                "",
                "- This is a bivariate robustness theorem for a candidate feature, not a strict dynamical source theorem.",
                "- No `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, L_total promotion, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
