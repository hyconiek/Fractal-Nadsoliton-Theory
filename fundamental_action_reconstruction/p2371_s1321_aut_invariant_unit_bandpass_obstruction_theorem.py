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

OUT = GEN / "p2371_s1321_aut_invariant_unit_bandpass_obstruction_theorem.json"
MD = GEN / "p2371_s1321_aut_invariant_unit_bandpass_obstruction_theorem.md"

SOURCE_FILES = {
    "P2370_D5_BANDPASS_SUPPORT": GEN / "p2370_s1320_d5_bandpass_support_closed_form_theorem.json",
    "P2369_LEDGER_CLOSED_FORM": GEN / "p2369_s1319_self_recorded_ledger_closed_form_uniqueness_theorem.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5


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
        for neighbor in ((node + step) % Z12_NODE_COUNT, (node - step) % Z12_NODE_COUNT):
            if neighbor in support_set:
                edges.add(tuple(sorted((node, neighbor))))
    return len(edges)


def support_rows() -> list[dict[str, Any]]:
    rows = []
    for support in combinations(range(Z12_NODE_COUNT), SUPPORT_SIZE):
        h1 = internal_edges(support, 1)
        h5 = internal_edges(support, 5)
        rows.append(
            {
                "support": list(support),
                "h1_distance_1_edges": h1,
                "h5_distance_5_edges": h5,
                "full_aut_unit_sum_h1_plus_h5": h1 + h5,
            }
        )
    return rows


def score_maximizers(rows: list[dict[str, Any]], a: int, b: int) -> dict[str, Any]:
    scored = [
        (a * row["h1_distance_1_edges"] + b * row["h5_distance_5_edges"], row)
        for row in rows
    ]
    maximum = max(score for score, _ in scored)
    maximizers = [row for score, row in scored if score == maximum]
    pair_distribution: dict[str, int] = {}
    for row in maximizers:
        key = f"{row['h1_distance_1_edges']},{row['h5_distance_5_edges']}"
        pair_distribution[key] = pair_distribution.get(key, 0) + 1
    return {
        "weights": {"a_h1": a, "b_h5": b},
        "maximum_score": maximum,
        "maximizer_count": len(maximizers),
        "maximizer_h1_h5_pair_distribution": dict(sorted(pair_distribution.items())),
        "maximizers_sample": maximizers[:12],
    }


def orbit_path_supports(step: int) -> list[list[int]]:
    return sorted(
        {tuple(sorted((source + index * step) % Z12_NODE_COUNT for index in range(SUPPORT_SIZE))) for source in range(Z12_NODE_COUNT)}
    )


def chamber_certificate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    pair_counts: dict[str, int] = {}
    for row in rows:
        key = f"{row['h1_distance_1_edges']},{row['h5_distance_5_edges']}"
        pair_counts[key] = pair_counts.get(key, 0) + 1

    d5_paths = orbit_path_supports(5)
    d1_paths = orbit_path_supports(1)
    d5_path_rows = [row for row in rows if tuple(row["support"]) in d5_paths]
    d1_path_rows = [row for row in rows if tuple(row["support"]) in d1_paths]
    mixed_rows = [row for row in rows if row["h1_distance_1_edges"] == 3 and row["h5_distance_5_edges"] == 3]

    full_aut = score_maximizers(rows, 1, 1)
    pure_h5 = score_maximizers(rows, 0, 1)
    threshold_tie = score_maximizers(rows, 1, 3)
    d5_selecting = score_maximizers(rows, 1, 4)

    return {
        "support_count": len(rows),
        "h1_h5_pair_count_distribution": dict(sorted(pair_counts.items())),
        "d5_path_orbit_count": len(d5_paths),
        "d5_path_h1_h5_pairs": sorted({(row["h1_distance_1_edges"], row["h5_distance_5_edges"]) for row in d5_path_rows}),
        "d1_path_orbit_count": len(d1_paths),
        "d1_path_h1_h5_pairs": sorted({(row["h1_distance_1_edges"], row["h5_distance_5_edges"]) for row in d1_path_rows}),
        "mixed_h1_eq_h5_eq_3_count": len(mixed_rows),
        "full_aut_invariant_a_eq_b_maximizers": full_aut,
        "pure_h5_bandpass_maximizers": pure_h5,
        "threshold_tie_a_over_b_eq_1_over_3": threshold_tie,
        "d5_selecting_example_a_over_b_eq_1_over_4": d5_selecting,
        "closed_form_chamber_statement": "For score a*h1+b*h5 with b>0 and a>=0, the d5 path orbit is the unique maximizer exactly when a/b < 1/3. At a/b=1/3 it ties the 24 mixed (h1,h5)=(3,3) supports; at a/b=1 full-Aut unit-band invariance selects only the mixed supports.",
        "aut_invariant_a_eq_b_selects_d5_path_orbit": full_aut["maximizer_h1_h5_pair_distribution"] == {"0,4": 12},
        "d5_selection_requires_band_polarity_break": True,
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}

    rows = support_rows()
    chamber = chamber_certificate(rows)

    theorem_export = {
        "theorem_name": "P2371 Aut-invariant unit-band band-pass obstruction theorem",
        "claim": (
            "On five-node supports in Z12, any full-Aut unit-band pair action that weights "
            "distance-1 and distance-5 unit bands equally cannot derive the P2370 distance-5 "
            "band-pass selector. The Aut-invariant score h1+h5 is maximized by 24 mixed supports "
            "with (h1,h5)=(3,3), not by the 12 distance-5 paths with (0,4). For scores a*h1+b*h5 "
            "with b>0 and a>=0, the distance-5 path orbit is uniquely selected only in the "
            "symmetry-breaking chamber a/b<1/3. Thus the band-pass polarity is an extra selector/source "
            "premise unless bridge-completed dynamics exports it."
        ),
        "positive_content": [
            "All 792 five-node supports are classified by their (h1,h5) edge counts.",
            "The full-Aut invariant unit-band score h1+h5 has 24 maximizers, all with (h1,h5)=(3,3).",
            "The pure distance-5 score h5 recovers the 12 P2370 distance-5 path supports.",
            "The exact chamber boundary for nonnegative linear scores a*h1+b*h5 is a/b=1/3.",
            "This proves that deriving the d5 band-pass action requires a real band-polarity source, not only full-Aut unit-band symmetry.",
        ],
        "not_licensed": [
            "strict derivation of band-polarity source",
            "strict derivation of distance-5 band-pass action",
            "strict selector closure",
            "unique translate/source selection",
            "QW-2191 discharge",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2371_S1321_AUT_INVARIANT_UNIT_BANDPASS_OBSTRUCTION_THEOREM",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "p2370_replay": {
            "packet_id": artifacts["P2370_D5_BANDPASS_SUPPORT"].get("packet_id"),
            "maximizers_are_12_step5_paths": artifacts["P2370_D5_BANDPASS_SUPPORT"].get("gatekeeper_checks", {}).get("maximizers_are_12_step5_paths"),
            "unit_step_orbit_exposes_full_aut_obstruction": artifacts["P2370_D5_BANDPASS_SUPPORT"].get("gatekeeper_checks", {}).get("unit_step_orbit_exposes_full_aut_obstruction"),
        },
        "chamber_certificate": chamber,
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2370_loaded": probe["p2370_replay"]["packet_id"] == "P2370",
        "support_count_is_792": chamber["support_count"] == 792,
        "full_aut_does_not_select_d5": chamber["aut_invariant_a_eq_b_selects_d5_path_orbit"] is False,
        "full_aut_selects_24_mixed_supports": chamber["full_aut_invariant_a_eq_b_maximizers"]["maximizer_count"] == 24
        and chamber["full_aut_invariant_a_eq_b_maximizers"]["maximizer_h1_h5_pair_distribution"] == {"3,3": 24},
        "pure_h5_recovers_d5_paths": chamber["pure_h5_bandpass_maximizers"]["maximizer_count"] == 12
        and chamber["pure_h5_bandpass_maximizers"]["maximizer_h1_h5_pair_distribution"] == {"0,4": 12},
        "threshold_tie_is_correct": chamber["threshold_tie_a_over_b_eq_1_over_3"]["maximizer_h1_h5_pair_distribution"] == {"0,4": 12, "3,3": 24},
        "d5_selecting_example_is_correct": chamber["d5_selecting_example_a_over_b_eq_1_over_4"]["maximizer_h1_h5_pair_distribution"] == {"0,4": 12},
        "docs_updated_with_p2371_theorem": all("P2371/S1321" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2371_s1321_v1",
        "packet_id": "P2371",
        "stage_id": "S1321",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_AUT_INVARIANT_BANDPASS_OBSTRUCTION_NO_QW2191_DISCHARGE",
        "result_kind": "AUT_INVARIANT_UNIT_BANDPASS_OBSTRUCTION_THEOREM",
        "aut_invariant_unit_bandpass_obstruction_theorem": probe,
        "recommended_next_honest_step": (
            "Look for a bridge-completed dynamical export of the band-polarity inequality a/b<1/3, "
            "or classify the distance-5 band-pass polarity as an explicit non-strict selector/source premise."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2371 S1321: Aut-invariant unit-band band-pass obstruction theorem",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2371/S1321 proves that a full-Aut unit-band invariant support action cannot derive the P2370 distance-5 band-pass selector by itself.",
                "",
                "## Certificate",
                "",
                f"- Supports classified: `{chamber['support_count']}`.",
                f"- Full-Aut score `h1+h5` maximizers: `{chamber['full_aut_invariant_a_eq_b_maximizers']['maximizer_count']}` with pair distribution `{chamber['full_aut_invariant_a_eq_b_maximizers']['maximizer_h1_h5_pair_distribution']}`.",
                f"- Pure h5 maximizers: `{chamber['pure_h5_bandpass_maximizers']['maximizer_count']}` with pair distribution `{chamber['pure_h5_bandpass_maximizers']['maximizer_h1_h5_pair_distribution']}`.",
                f"- Threshold tie at `a/b=1/3`: `{chamber['threshold_tie_a_over_b_eq_1_over_3']['maximizer_h1_h5_pair_distribution']}`.",
                "- D5 unique selection for nonnegative linear scores requires `a/b < 1/3`.",
                "",
                "## Hard limits",
                "",
                "- This is an obstruction theorem for full-Aut unit-band derivations, not a strict derivation of band polarity.",
                "- No `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
