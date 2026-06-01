#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from itertools import combinations
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"

OUT = GEN / "p2385_s1335_exact_z12_support_chamber_theorem.json"
MD = GEN / "p2385_s1335_exact_z12_support_chamber_theorem.md"

SOURCE_FILES = {
    "P2384_SYMBOLIC_INEQUALITY": GEN / "p2384_s1334_symbolic_bathtub_inequality_proof_packet.json",
    "P2382_BOUNDED_DENSITY_BATHTUB": GEN / "p2382_s1332_bounded_density_bathtub_frontload_certificate.json",
    "P2370_D5_BANDPASS_SUPPORT": GEN / "p2370_s1320_d5_bandpass_support_closed_form_theorem.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
D5_STEP = 5
H1_STEP = 1
CHAMBER_RATIO_NUMERATOR = 1
CHAMBER_RATIO_DENOMINATOR = 3


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


def internal_edges(support: tuple[int, ...], step: int) -> int:
    support_set = set(support)
    edges = set()
    for node in support_set:
        neighbor = (node + step) % Z12_NODE_COUNT
        if neighbor in support_set:
            edges.add(tuple(sorted((node, neighbor))))
    return len(edges)


def support_rows() -> list[dict[str, Any]]:
    rows = []
    for support in combinations(range(Z12_NODE_COUNT), SUPPORT_SIZE):
        h1 = internal_edges(support, H1_STEP)
        h5 = internal_edges(support, D5_STEP)
        rows.append({"support": list(support), "h1": h1, "h5": h5})
    return rows


def d5_cycle_path_supports() -> list[list[int]]:
    supports = []
    for start in range(Z12_NODE_COUNT):
        support = sorted((start + i * D5_STEP) % Z12_NODE_COUNT for i in range(SUPPORT_SIZE))
        supports.append(support)
    return sorted(supports)


def pair_distribution(rows: list[dict[str, Any]]) -> dict[str, int]:
    distribution: dict[str, int] = {}
    for row in rows:
        key = f"{row['h1']},{row['h5']}"
        distribution[key] = distribution.get(key, 0) + 1
    return dict(sorted(distribution.items(), key=lambda item: tuple(map(int, item[0].split(",")))))


def chamber_gap_table(distribution: dict[str, int]) -> list[dict[str, Any]]:
    table = []
    for key, count in distribution.items():
        h1, h5 = (int(part) for part in key.split(","))
        numerator_at_boundary = CHAMBER_RATIO_DENOMINATOR * (4 - h5) - h1
        target_pair = h1 == 0 and h5 == 4
        table.append(
            {
                "h1": h1,
                "h5": h5,
                "count": count,
                "target_pair_0_4": target_pair,
                "score_gap_from_0_4_over_b_at_r_equals_1_3_numerator": numerator_at_boundary,
                "strictly_beaten_for_r_lt_1_3": target_pair or numerator_at_boundary >= 0,
                "boundary_tie_pair_at_r_1_3": (not target_pair) and numerator_at_boundary == 0,
            }
        )
    return table


def exact_support_chamber(rows: list[dict[str, Any]]) -> dict[str, Any]:
    distribution = pair_distribution(rows)
    d5_paths = d5_cycle_path_supports()
    enumerated_targets = sorted(row["support"] for row in rows if row["h1"] == 0 and row["h5"] == 4)
    table = chamber_gap_table(distribution)
    boundary_ties = [row for row in table if row["boundary_tie_pair_at_r_1_3"]]
    invalid_gap_rows = [row for row in table if not row["strictly_beaten_for_r_lt_1_3"]]
    max_h5 = max(row["h5"] for row in rows)
    max_h5_rows = [row for row in rows if row["h5"] == max_h5]
    max_h5_distribution = pair_distribution(max_h5_rows)
    return {
        "total_supports_checked": len(rows),
        "expected_support_count_binomial_12_5": 792,
        "pair_distribution": distribution,
        "max_h5": max_h5,
        "max_h5_pair_distribution": max_h5_distribution,
        "d5_path_supports": d5_paths,
        "enumerated_target_supports": enumerated_targets,
        "d5_path_supports_match_enumerated_targets": d5_paths == enumerated_targets,
        "chamber_ratio_condition": "b>0, a>=0, a/b<1/3",
        "proof_rule": (
            "The (h1,h5)=(0,4) supports have score 4b. Every other pair has score gap "
            "4b-(a*h1+b*h5)=b*((4-h5)-r*h1), r=a/b. At r=1/3 the integer numerator "
            "3*(4-h5)-h1 is nonnegative for every observed pair; the only non-target zero pair is (3,3), "
            "so the inequality is strict for r<1/3."
        ),
        "gap_table": table,
        "boundary_ties_at_r_1_3": boundary_ties,
        "invalid_gap_rows": invalid_gap_rows,
        "target_count": len(enumerated_targets),
        "target_pair_distribution": {"0,4": len(enumerated_targets)},
    }


def p2382_replay_summary(artifact: dict[str, Any]) -> dict[str, Any]:
    try:
        cert = artifact["bounded_density_bathtub_frontload_theorem"]["bounded_density_bathtub_frontload_certificate"]
        obligation = cert["rectangle_worst_cap_source_obligation"]
        replay = obligation["cap_test_support_score"]
        below = obligation["below_threshold_support_score"]
        return {
            "cap_threshold_gt": obligation["cap_threshold_gt"],
            "cap_test_a_over_b": replay["a_over_b"],
            "cap_test_d5_chamber": replay["d5_chamber"],
            "cap_test_pair_distribution": replay["maximizer_h1_h5_pair_distribution"],
            "below_threshold_a_over_b": below["a_over_b"],
            "below_threshold_d5_chamber": below["d5_chamber"],
        }
    except KeyError:
        return {"status": "P2382_REPLAY_FIELDS_MISSING"}


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}

    rows = support_rows()
    chamber = exact_support_chamber(rows)
    p2382_summary = p2382_replay_summary(artifacts["P2382_BOUNDED_DENSITY_BATHTUB"])

    theorem_export = {
        "name": "P2385/S1335 exact Z12 support chamber theorem",
        "claim": (
            "For 5-subsets of Z12 scored by a*h1+b*h5, the chamber b>0, a>=0, a/b<1/3 has exactly the 12 step-5 path supports as unique maximizers. "
            "Thus the P2382/P2384 analytic cap proof can hand off to an exact finite support theorem rather than to an unstructured enumeration."
        ),
        "positive_content": [
            "Classifies the full (h1,h5) pair distribution for all binomial(12,5)=792 supports.",
            "Exports an integer gap table proving the chamber inequality against the (0,4) target pair.",
            "Identifies the 12 target supports as length-5 paths in the Z12 step-5 cycle.",
        ],
        "not_licensed": [
            "strict source theorem deriving the cap M or bang-bang density",
            "promotion of rho_M, W_M, or support selection into L_total",
            "claim that the cap/frontload premise is strict-core supplied",
            "QW-2191 discharge or selector closure",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2385_S1335_EXACT_Z12_SUPPORT_CHAMBER_THEOREM",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "artifact_replay": {
            "p2384_packet_id": artifacts["P2384_SYMBOLIC_INEQUALITY"].get("packet_id"),
            "p2382_packet_id": artifacts["P2382_BOUNDED_DENSITY_BATHTUB"].get("packet_id"),
            "p2370_packet_id": artifacts["P2370_D5_BANDPASS_SUPPORT"].get("packet_id"),
        },
        "exact_support_chamber_certificate": chamber,
        "p2382_cap_chamber_replay_summary": p2382_summary,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2384_loaded": probe["artifact_replay"]["p2384_packet_id"] == "P2384",
        "p2382_loaded": probe["artifact_replay"]["p2382_packet_id"] == "P2382",
        "p2370_loaded": probe["artifact_replay"]["p2370_packet_id"] == "P2370",
        "all_792_supports_checked": chamber["total_supports_checked"] == chamber["expected_support_count_binomial_12_5"],
        "target_supports_are_exact_d5_paths": chamber["d5_path_supports_match_enumerated_targets"],
        "target_count_is_12": chamber["target_count"] == 12,
        "max_h5_pair_is_only_0_4": chamber["max_h5_pair_distribution"] == {"0,4": 12},
        "no_invalid_chamber_gap_rows": chamber["invalid_gap_rows"] == [],
        "only_boundary_tie_pair_is_3_3": chamber["boundary_ties_at_r_1_3"] == [
            {
                "h1": 3,
                "h5": 3,
                "count": 24,
                "target_pair_0_4": False,
                "score_gap_from_0_4_over_b_at_r_equals_1_3_numerator": 0,
                "strictly_beaten_for_r_lt_1_3": True,
                "boundary_tie_pair_at_r_1_3": True,
            }
        ],
        "p2382_cap_test_enters_exact_chamber": p2382_summary.get("cap_test_d5_chamber") is True
        and p2382_summary.get("cap_test_pair_distribution") == {"0,4": 12},
        "p2382_negative_control_still_fails_chamber": p2382_summary.get("below_threshold_d5_chamber") is False,
        "docs_updated_with_p2385_support_chamber": all("P2385/S1335" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2385_s1335_v1",
        "packet_id": "P2385",
        "stage_id": "S1335",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_EXACT_Z12_SUPPORT_CHAMBER_THEOREM_SOURCE_STILL_OPEN",
        "result_kind": "EXACT_Z12_SUPPORT_CHAMBER_THEOREM",
        "exact_z12_support_chamber_theorem": probe,
        "recommended_next_honest_step": (
            "The analytic cap proof and the exact finite support chamber are now separated. The remaining honest bottleneck is a real strict source theorem for the bounded front-loaded density/cap, not more support enumeration."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2385 S1335: exact Z12 support chamber theorem",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2385/S1335 separates the finite support-selection theorem from the analytic cap proof.",
                "For 5-subsets of `Z12` scored by `a*h1+b*h5`, the chamber",
                "",
                "```text",
                "b>0, a>=0, a/b<1/3",
                "```",
                "",
                "has exactly the 12 `(h1,h5)=(0,4)` supports as unique maximizers.",
                "These supports are precisely the length-5 paths in the step-5 cycle on `Z12`.",
                "",
                "## Finite certificate",
                "",
                f"- Supports checked: `{chamber['total_supports_checked']}`.",
                f"- Full pair distribution: `{chamber['pair_distribution']}`.",
                f"- Target support count: `{chamber['target_count']}`.",
                f"- Boundary tie pair at `a/b=1/3`: `{chamber['boundary_ties_at_r_1_3']}`.",
                "",
                "## Hard limits",
                "",
                "- This is an exact finite chamber theorem for the support-scoring layer, not a strict source theorem for the cap or bang-bang profile.",
                "- No `L_total` promotion, `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
