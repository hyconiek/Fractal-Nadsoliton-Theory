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
SCRATCH = ROOT / "scratch"

OUT = GEN / "p2370_s1320_d5_bandpass_support_closed_form_theorem.json"
MD = GEN / "p2370_s1320_d5_bandpass_support_closed_form_theorem.md"

SOURCE_FILES = {
    "P2369_LEDGER_CLOSED_FORM": GEN
    / "p2369_s1319_self_recorded_ledger_closed_form_uniqueness_theorem.json",
    "P2368_ENDPOINT_ANCHOR": GEN
    / "p2368_s1318_self_recorded_endpoint_anchor_selector_candidate_probe.json",
    "SCRATCH_BANDPASS_D5": SCRATCH
    / "bridge_strict_alpha_bandpass_d5_support_selector_certificate_report.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
BANDPASS_STEP = 5
BALANCED_LEDGER = (2, 2, 2, 1, 1)


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


def step_neighbors(node: int, step: int = BANDPASS_STEP) -> tuple[int, int]:
    return ((node + step) % Z12_NODE_COUNT, (node - step) % Z12_NODE_COUNT)


def ordered_step_path(source: int, step: int = BANDPASS_STEP, size: int = SUPPORT_SIZE) -> tuple[int, ...]:
    return tuple((source + index * step) % Z12_NODE_COUNT for index in range(size))


def canonical_support(nodes: tuple[int, ...] | set[int]) -> tuple[int, ...]:
    return tuple(sorted(nodes))


def internal_step_edges(support: tuple[int, ...], step: int = BANDPASS_STEP) -> list[tuple[int, int]]:
    support_set = set(support)
    edges = set()
    for node in support_set:
        for neighbor in step_neighbors(node, step):
            if neighbor in support_set:
                edges.add(tuple(sorted((node, neighbor))))
    return sorted(edges)


def connected_components(support: tuple[int, ...], step: int = BANDPASS_STEP) -> list[list[int]]:
    support_set = set(support)
    seen: set[int] = set()
    components: list[list[int]] = []
    for start in sorted(support_set):
        if start in seen:
            continue
        stack = [start]
        component = []
        seen.add(start)
        while stack:
            node = stack.pop()
            component.append(node)
            for neighbor in step_neighbors(node, step):
                if neighbor in support_set and neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        components.append(sorted(component))
    return components


def cyclic_gap_profile_in_step_order(support: tuple[int, ...], step: int = BANDPASS_STEP) -> tuple[int, ...]:
    inverse = pow(step, -1, Z12_NODE_COUNT)
    coordinates = sorted((node * inverse) % Z12_NODE_COUNT for node in support)
    gaps = []
    for left, right in zip(coordinates, coordinates[1:] + [coordinates[0] + Z12_NODE_COUNT]):
        gaps.append(right - left)
    return tuple(sorted(gaps))


def bandpass_closed_form_certificate() -> dict[str, Any]:
    all_supports = [tuple(row) for row in combinations(range(Z12_NODE_COUNT), SUPPORT_SIZE)]
    rows = []
    h_counts: dict[int, int] = {}
    for support in all_supports:
        edges = internal_step_edges(support)
        h5 = len(edges)
        h_counts[h5] = h_counts.get(h5, 0) + 1
        components = connected_components(support)
        rows.append(
            {
                "support": list(support),
                "h5_edges": h5,
                "component_count": len(components),
                "component_sizes": [len(component) for component in components],
                "cyclic_gap_profile_in_step_order": list(cyclic_gap_profile_in_step_order(support)),
            }
        )
    max_h5 = max(h_counts)
    maximizers = [row for row in rows if row["h5_edges"] == max_h5]
    expected_path_supports = sorted({canonical_support(ordered_step_path(source)) for source in range(Z12_NODE_COUNT)})
    maximizer_supports = sorted(tuple(row["support"]) for row in maximizers)
    return {
        "cycle_statement": "Since gcd(5,12)=1, distance-5 edges form a single 12-cycle. Any 5-vertex induced subgraph of a cycle has at most 4 internal edges; equality holds exactly when the selected vertices are one connected 5-path in that cycle.",
        "gcd_5_12": 1,
        "support_count_cross_check": len(all_supports),
        "h5_count_distribution": {str(key): h_counts[key] for key in sorted(h_counts)},
        "closed_form_max_h5": SUPPORT_SIZE - 1,
        "brute_force_max_h5": max_h5,
        "maximizer_count": len(maximizers),
        "expected_path_support_count": len(expected_path_supports),
        "maximizers_equal_step5_path_orbit": maximizer_supports == expected_path_supports,
        "maximizers": maximizers,
        "expected_path_supports": [list(row) for row in expected_path_supports],
        "all_maximizers_connected_5_paths": all(
            row["component_count"] == 1 and row["component_sizes"] == [SUPPORT_SIZE]
            for row in maximizers
        ),
    }


def unit_step_orbit_certificate() -> dict[str, Any]:
    units = [1, 5, 7, 11]
    images_of_step5 = sorted({(unit * BANDPASS_STEP) % Z12_NODE_COUNT for unit in units})
    unoriented_images = sorted({min(step, (-step) % Z12_NODE_COUNT) for step in images_of_step5})
    return {
        "aut_z12_units": units,
        "oriented_images_of_step5": images_of_step5,
        "unoriented_images_of_step5": unoriented_images,
        "step5_unoriented_pair_fixed_by_aut_z12": unoriented_images == [5],
        "step5_unoriented_pair_collapses_with_step1_under_full_aut_z12": unoriented_images == [1, 5],
        "interpretation": "Full Aut(Z12) sends the unoriented distance-5 band into the unit-band pair {1,5}; therefore the distance-5 band-pass action is itself an extra support/source premise, not a full-Aut strict-core consequence.",
    }


def conditional_chain_certificate(bandpass: dict[str, Any], p2369: dict[str, Any]) -> dict[str, Any]:
    p2369_probe = p2369.get("self_recorded_ledger_closed_form_uniqueness_theorem", {})
    return {
        "premise_1_bandpass_support": "Admit a distance-5 band-pass support action maximizing h5 internal edges.",
        "consequence_1_support_orbit": "The maximizers are exactly the 12 translate supports of a 5-node path in the distance-5 cycle.",
        "premise_2_ordered_path_traversal": "Choose one of the two orientations of that path only as an endpoint-ordering convention, not as strict selector closure.",
        "consequence_2_ledger": "P2369 closed-form proof uniquely selects the ordered ledger (2,2,2,1,1) by lexicographic ripple/arrow.",
        "consequence_3_endpoint_anchor": "Given the ordered d5 path and ledger, P2368 recovers source/orientation equivariantly from endpoint values.",
        "p2369_loaded_packet": p2369.get("packet_id"),
        "p2369_arrow_tiebreak_unique": p2369.get("gatekeeper_checks", {}).get("arrow_tiebreak_unique"),
        "p2369_endpoint_anchor_distinct": p2369.get("gatekeeper_checks", {}).get("endpoint_anchor_distinct"),
        "bandpass_maximizers_equal_step5_orbit": bandpass["maximizers_equal_step5_path_orbit"],
        "remaining_missing_premise": "strict derivation of the distance-5 band-pass action and the path-order/arrow-action convention from bridge-completed nadsoliton dynamics",
        "p2369_theorem_fingerprint": p2369_probe.get("theorem_fingerprint_sha256"),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}

    bandpass = bandpass_closed_form_certificate()
    unit_orbit = unit_step_orbit_certificate()
    chain = conditional_chain_certificate(bandpass, artifacts["P2369_LEDGER_CLOSED_FORM"])

    theorem_export = {
        "theorem_name": "P2370 conditional distance-5 band-pass support theorem",
        "claim": (
            "If a distance-5 band-pass support action is admitted, then the support side of the "
            "P2368/P2369 selector candidate has a closed-form finite proof: the distance-5 graph "
            "is a 12-cycle, every five-vertex induced subgraph has at most four distance-5 edges, "
            "and equality occurs exactly on the twelve translated five-node paths. Full Aut(Z12) still mixes "
            "the distance-5 band with the distance-1 unit band, so combined with P2369 this leaves "
            "the strict derivation of the band-pass action and arrow/path-order "
            "convention as the real remaining selector premise."
        ),
        "positive_content": [
            "Closed-form cycle argument proves max h5=4 for five-node supports.",
            "A 792-support brute-force scan is retained only as a cross-check and agrees with the closed-form theorem.",
            "The maximizers are exactly the 12 translated distance-5 path supports.",
            "Full Aut(Z12) does not preserve the distance-5 band alone; it mixes the unoriented unit bands {1,5}, so the band-pass action remains a visible extra premise.",
            "Together with P2369, finite support and finite ledger uniqueness are isolated from the still-open strict source premise.",
        ],
        "not_licensed": [
            "strict derivation of the distance-5 band-pass action",
            "strict derivation of path orientation/order convention",
            "strict derivation of the arrow action from nadsoliton dynamics",
            "unique translate/source selection",
            "QW-2191 discharge",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2370_S1320_D5_BANDPASS_SUPPORT_CLOSED_FORM_THEOREM",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "bandpass_closed_form_certificate": bandpass,
        "unit_step_orbit_certificate": unit_orbit,
        "conditional_chain_certificate": chain,
        "scratch_bandpass_replay": {
            "result_kind": artifacts["SCRATCH_BANDPASS_D5"].get("result_kind"),
            "max_h5": artifacts["SCRATCH_BANDPASS_D5"].get("bandpass_selector_scan", {}).get("max_h5"),
            "maximizer_count": artifacts["SCRATCH_BANDPASS_D5"].get("bandpass_selector_scan", {}).get("maximizer_count"),
            "maximizers_equal_fifth_step_orbit": artifacts["SCRATCH_BANDPASS_D5"].get("bandpass_selector_scan", {}).get("maximizers_equal_fifth_step_orbit"),
        },
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2369_loaded": artifacts["P2369_LEDGER_CLOSED_FORM"].get("packet_id") == "P2369",
        "closed_form_max_matches_bruteforce": bandpass["closed_form_max_h5"] == bandpass["brute_force_max_h5"],
        "support_count_is_792": bandpass["support_count_cross_check"] == 792,
        "maximizers_are_12_step5_paths": bandpass["maximizer_count"] == 12 and bandpass["maximizers_equal_step5_path_orbit"],
        "all_maximizers_connected_paths": bandpass["all_maximizers_connected_5_paths"],
        "unit_step_orbit_exposes_full_aut_obstruction": unit_orbit["step5_unoriented_pair_collapses_with_step1_under_full_aut_z12"]
        and not unit_orbit["step5_unoriented_pair_fixed_by_aut_z12"],
        "scratch_bandpass_replay_agrees": probe["scratch_bandpass_replay"] == {
            "result_kind": "SCRATCH_STRICT_ALPHA_BANDPASS_D5_SUPPORT_SELECTOR_CERTIFICATE_PROBE__NOT_A_THEOREM",
            "max_h5": 4,
            "maximizer_count": 12,
            "maximizers_equal_fifth_step_orbit": True,
        },
        "docs_updated_with_p2370_theorem": all("P2370/S1320" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2370_s1320_v1",
        "packet_id": "P2370",
        "stage_id": "S1320",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_D5_BANDPASS_SUPPORT_CLOSED_FORM_NO_QW2191_DISCHARGE",
        "result_kind": "D5_BANDPASS_SUPPORT_CLOSED_FORM_THEOREM",
        "d5_bandpass_support_closed_form_theorem": probe,
        "recommended_next_honest_step": (
            "Try to derive the distance-5 band-pass action itself from bridge-completed nadsoliton dynamics "
            "or prove that the band-pass action is an extra non-strict selector/source premise; do not treat "
            "the finite support theorem as QW-2191 discharge."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2370 S1320: distance-5 band-pass support closed-form theorem",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2370/S1320 proves the conditional support side of the self-recorded selector chain: if a distance-5 band-pass support action is admitted, its five-node maximizers are exactly the twelve translated paths in the distance-5 cycle.",
                "",
                "## Certificate",
                "",
                f"- Supports checked as cross-check: `{bandpass['support_count_cross_check']}`.",
                f"- Closed-form max h5: `{bandpass['closed_form_max_h5']}`; brute-force max h5: `{bandpass['brute_force_max_h5']}`.",
                f"- Maximizer count: `{bandpass['maximizer_count']}`.",
                f"- Maximizers equal step-5 path orbit: `{bandpass['maximizers_equal_step5_path_orbit']}`.",
                f"- All maximizers connected 5-paths: `{bandpass['all_maximizers_connected_5_paths']}`.",
                f"- Full Aut(Z12) mixes distance-5 with distance-1 unit bands: `{unit_orbit['step5_unoriented_pair_collapses_with_step1_under_full_aut_z12']}`.",
                "",
                "## Hard limits",
                "",
                "- This proves only the conditional finite support theorem after admitting a distance-5 band-pass action.",
                "- It does not derive the band-pass action, path-order convention, arrow action, unique translate/source, QW-2191 discharge, legacy role transfer, or ToE closure.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
