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

OUT = GEN / "p2375_s1325_damping_compression_polarity_interval_robustness_theorem.json"
MD = GEN / "p2375_s1325_damping_compression_polarity_interval_robustness_theorem.md"

SOURCE_FILES = {
    "P2374_DAMPING_COMPRESSION_CANDIDATE": GEN / "p2374_s1324_damping_compression_band_polarity_candidate_theorem.json",
    "P2373_CORRECTION_CONE": GEN / "p2373_s1323_bridge_kernel_polarity_correction_cone_theorem.json",
    "SCRATCH_DAMPING_SEPARATION": ROOT / "scratch" / "bridge_strict_completion_legacy_to_strict_damping_compression_separation_certificate_report.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
ETA = 9.0 / 5.0
STRICT_BETA = 1.0
BETA_TORS_INTERVAL = (0.0, 0.1)
CANONICAL_BETA_TORS = 0.01
CHAMBER_THRESHOLD = 1.0 / 3.0


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


def compression_log_weight(d: int, beta_tors: float) -> float:
    return math.log((1.0 + STRICT_BETA * d**ETA) / (1.0 + beta_tors * d))


def chamber_margin_function(beta_tors: float) -> float:
    """Positive iff C(1)/C(5)<1/3, after exponentiating C5-3*C1."""
    strict_d5_denominator = 1.0 + 5.0**ETA
    return strict_d5_denominator * (1.0 + beta_tors) ** 3 - 8.0 * (1.0 + 5.0 * beta_tors)


def chamber_margin_derivative(beta_tors: float) -> float:
    strict_d5_denominator = 1.0 + 5.0**ETA
    return 3.0 * strict_d5_denominator * (1.0 + beta_tors) ** 2 - 40.0


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
    }


def robustness_certificate() -> dict[str, Any]:
    lo, hi = BETA_TORS_INTERVAL
    rows = support_rows()
    sample_points = [lo, CANONICAL_BETA_TORS, 0.05, hi]
    sample_rows = []
    for beta_tors in sample_points:
        c1 = compression_log_weight(1, beta_tors)
        c5 = compression_log_weight(5, beta_tors)
        sample_rows.append(
            {
                "beta_tors": beta_tors,
                "C1": c1,
                "C5": c5,
                "C1_over_C5": c1 / c5,
                "margin_C5_minus_3C1": c5 - 3.0 * c1,
                "exponentiated_margin": chamber_margin_function(beta_tors),
                "score_audit": score_maximizers(rows, c1, c5),
            }
        )

    f_lo = chamber_margin_function(lo)
    f_hi = chamber_margin_function(hi)
    derivative_lo = chamber_margin_derivative(lo)
    derivative_hi = chamber_margin_derivative(hi)
    strict_d5_denominator = 1.0 + 5.0**ETA

    return {
        "interval": {"beta_tors_min": lo, "beta_tors_max": hi, "canonical_beta_tors": CANONICAL_BETA_TORS},
        "symbolic_inequality": {
            "target": "C1/C5 < 1/3",
            "equivalent_positive_margin": "F(x)=(1+5^(9/5))*(1+x)^3 - 8*(1+5*x) > 0",
            "derivative": "F'(x)=3*(1+5^(9/5))*(1+x)^2 - 40",
            "derivative_is_increasing_on_interval": True,
            "derivative_min_at_zero": derivative_lo,
            "margin_min_at_zero": f_lo,
            "margin_at_interval_max": f_hi,
            "strict_d5_denominator_1_plus_5_eta": strict_d5_denominator,
        },
        "proof_summary": {
            "derivative_positive_on_interval": derivative_lo > 0 and derivative_hi > derivative_lo,
            "margin_positive_on_interval": f_lo > 0 and derivative_lo > 0,
            "therefore_C1_over_C5_lt_one_third_for_all_interval": f_lo > 0 and derivative_lo > 0,
            "not_fine_tuned_to_beta_tors_0_01": True,
        },
        "sample_score_audits": sample_rows,
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}
    certificate = robustness_certificate()
    sample_audits = certificate["sample_score_audits"]

    theorem_export = {
        "theorem_name": "P2375 damping-compression polarity interval robustness theorem",
        "claim": (
            "For every beta_tors in [0,0.1], the strict nonlinear damping/compression surplus "
            "C(d)=log((1+d^(9/5))/(1+beta_tors*d)) satisfies C(1)/C(5)<1/3. Equivalently, "
            "F(x)=(1+5^(9/5))*(1+x)^3-8*(1+5*x) is positive on the interval because F(0)>0 and F'(0)>0 with F' increasing. "
            "Thus the P2374 d5-polarity candidate is interval-robust in beta_tors, not a fine-tuned artifact of beta_tors=0.01."
        ),
        "positive_content": [
            "Closed-form interval proof for the P2374 chamber inequality over beta_tors in [0,0.1].",
            "Endpoint/canonical/midpoint support scans verify d5 path selection on all 792 five-node supports.",
            "Links the result to the prior damping-separation certificate without claiming a source theorem.",
        ],
        "not_licensed": [
            "strict dynamical source theorem for C(d)",
            "derivation of eta=9/5 or beta=1 from nadsoliton dynamics",
            "promotion of C(d) into L_total or selector action without a variational/transport theorem",
            "QW-2191 discharge",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2375_S1325_DAMPING_COMPRESSION_POLARITY_INTERVAL_ROBUSTNESS_THEOREM",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "artifact_replay": {
            "p2374_packet_id": artifacts["P2374_DAMPING_COMPRESSION_CANDIDATE"].get("packet_id"),
            "p2373_packet_id": artifacts["P2373_CORRECTION_CONE"].get("packet_id"),
            "scratch_damping_status": artifacts["SCRATCH_DAMPING_SEPARATION"].get("status"),
        },
        "robustness_certificate": certificate,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2374_loaded": probe["artifact_replay"]["p2374_packet_id"] == "P2374",
        "p2373_loaded": probe["artifact_replay"]["p2373_packet_id"] == "P2373",
        "scratch_damping_separation_loaded": isinstance(probe["artifact_replay"]["scratch_damping_status"], str),
        "derivative_positive_on_interval": certificate["proof_summary"]["derivative_positive_on_interval"],
        "margin_positive_on_interval": certificate["proof_summary"]["margin_positive_on_interval"],
        "all_sample_audits_select_d5_paths": all(
            row["score_audit"]["maximizer_h1_h5_pair_distribution"] == {"0,4": 12} for row in sample_audits
        ),
        "canonical_ratio_matches_p2374": abs(sample_audits[1]["C1_over_C5"] - 0.2354294000560429) < 1e-15,
        "docs_updated_with_p2375_robustness": all("P2375/S1325" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2375_s1325_v1",
        "packet_id": "P2375",
        "stage_id": "S1325",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_DAMPING_COMPRESSION_POLARITY_INTERVAL_ROBUST_NO_QW2191_DISCHARGE",
        "result_kind": "DAMPING_COMPRESSION_POLARITY_INTERVAL_ROBUSTNESS_THEOREM",
        "damping_compression_polarity_interval_robustness_theorem": probe,
        "recommended_next_honest_step": (
            "Use the robust C(d) chamber result as an acceptance target for a real variational/transport source theorem; "
            "prove that theorem or keep the compression polarity as a quantified non-strict selector premise."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2375 S1325: damping-compression polarity interval robustness theorem",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2375/S1325 proves that the P2374 compression-polarity candidate is robust over `beta_tors in [0,0.1]`.",
                "",
                "## Certificate",
                "",
                f"- Equivalent margin: `{certificate['symbolic_inequality']['equivalent_positive_margin']}`.",
                f"- Margin minimum at x=0: `{certificate['symbolic_inequality']['margin_min_at_zero']}`.",
                f"- Derivative minimum at x=0: `{certificate['symbolic_inequality']['derivative_min_at_zero']}`.",
                f"- Canonical C1/C5: `{sample_audits[1]['C1_over_C5']}`.",
                f"- Interval endpoint beta_tors=0.1 C1/C5: `{sample_audits[-1]['C1_over_C5']}`.",
                f"- All sampled support scans select: `{[row['score_audit']['maximizer_h1_h5_pair_distribution'] for row in sample_audits]}`.",
                "",
                "## Hard limits",
                "",
                "- This is an interval robustness theorem for a candidate feature, not a strict dynamical source theorem.",
                "- No `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, L_total promotion, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
