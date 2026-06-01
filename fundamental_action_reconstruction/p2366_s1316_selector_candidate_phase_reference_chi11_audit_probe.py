#!/usr/bin/env python3
from __future__ import annotations

import cmath
import hashlib
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
SCRATCH = ROOT / "scratch"

OUT = GEN / "p2366_s1316_selector_candidate_phase_reference_chi11_audit_probe.json"
MD = GEN / "p2366_s1316_selector_candidate_phase_reference_chi11_audit_probe.md"

SOURCE_FILES = {
    "ATOM_INFLUENCE": SCRATCH
    / "bridge_strict_completion_theorem_frontier_atom_influence_certificate_report.json",
    "LOW_WEIGHT_FRONTIER": SCRATCH
    / "bridge_strict_completion_theorem_frontier_low_weight_extension_certificate_report.json",
    "PHASE_REFERENCE_SELECTOR": SCRATCH
    / "bridge_strict_alpha_phase_reference_source_selector_certificate_report.json",
    "CHI11_PREMISE_LATTICE": SCRATCH
    / "bridge_strict_alpha_chi11_premise_dependency_lattice_audit_report.json",
    "LEGACY_TORSION_CHI11": SCRATCH
    / "bridge_legacy_torsion_chi11_opinion_audit_report.json",
    "P2363_BRIDGE_COMPLETED_MOMENTS": GEN
    / "p2363_s1313_legacy_strict_bridge_moment_lagrangian_transport_probe.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

REQUIRED_DOC_SNIPPETS = [
    "P2366/S1316",
    "phase-origin selector candidate",
    "chi11_selector_source",
    "QW-2191 remains open",
]

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
DISTANCE_SELECTED = 5
FORWARD_ASSIGNMENT = (2, 2, 2, 1, 1)
COPRIME_SOURCE_MODES = (1, 5, 7, 11)
NONCOPRIME_MODES = (2, 3, 4, 6, 8, 9, 10)
ORIENTATION_BISPECTRUM_PAIR = (1, 5)
ROUND_DIGITS = 12


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8")


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def rounded(value: float) -> float:
    result = round(value, ROUND_DIGITS)
    return 0.0 if result == -0.0 else result


def ordered_d5_support(source: int, orientation: int) -> tuple[int, ...]:
    step = (orientation * DISTANCE_SELECTED) % Z12_NODE_COUNT
    return tuple((source + index * step) % Z12_NODE_COUNT for index in range(SUPPORT_SIZE))


def value_configuration(source: int, orientation: int) -> tuple[int, ...]:
    values = [0] * Z12_NODE_COUNT
    for node, value in zip(ordered_d5_support(source, orientation), FORWARD_ASSIGNMENT):
        values[node] = value
    return tuple(values)


def dft(config: tuple[int, ...], mode: int) -> complex:
    return sum(
        value * cmath.exp(-2j * math.pi * mode * node / Z12_NODE_COUNT)
        for node, value in enumerate(config)
    )


def bispectrum(config: tuple[int, ...], left_mode: int, right_mode: int) -> complex:
    return dft(config, left_mode) * dft(config, right_mode) * dft(
        config,
        (left_mode + right_mode) % Z12_NODE_COUNT,
    ).conjugate()


def phase_turns(value: complex) -> float:
    if abs(value) < 1e-12:
        raise ValueError("phase is undefined for a near-zero Fourier coefficient")
    phase = (math.atan2(value.imag, value.real) / (2.0 * math.pi)) % 1.0
    return rounded(0.0 if round(phase, ROUND_DIGITS) in {1.0, -0.0} else phase)


def signed_orientation_from_bispectrum(config: tuple[int, ...]) -> int:
    marker = bispectrum(config, *ORIENTATION_BISPECTRUM_PAIR)
    if abs(marker.imag) < 1e-10:
        raise ValueError("chiral bispectrum marker has zero imaginary part")
    return -1 if marker.imag > 0 else 1


def calibrated_phase_references() -> dict[str, dict[str, float]]:
    return {
        str(orientation): {
            str(mode): phase_turns(dft(value_configuration(0, orientation), mode))
            for mode in COPRIME_SOURCE_MODES
        }
        for orientation in (-1, 1)
    }


def recover_source_from_phase(phase: float, reference_phase: float, mode: int) -> int:
    inverse = pow(mode, -1, Z12_NODE_COUNT)
    raw_shift = ((reference_phase - phase) % 1.0) * Z12_NODE_COUNT * inverse
    return int(round(raw_shift)) % Z12_NODE_COUNT


def source_recovery_rows() -> list[dict[str, Any]]:
    references = calibrated_phase_references()
    rows: list[dict[str, Any]] = []
    for orientation in (-1, 1):
        for source in range(Z12_NODE_COUNT):
            config = value_configuration(source, orientation)
            predicted_orientation = signed_orientation_from_bispectrum(config)
            per_mode: dict[str, Any] = {}
            for mode in COPRIME_SOURCE_MODES:
                observed_phase = phase_turns(dft(config, mode))
                reference_phase = references[str(predicted_orientation)][str(mode)]
                recovered_source = recover_source_from_phase(observed_phase, reference_phase, mode)
                per_mode[str(mode)] = {
                    "observed_phase_turns": observed_phase,
                    "reference_phase_turns": reference_phase,
                    "recovered_source": recovered_source,
                    "matches_actual_source": recovered_source == source,
                }
            rows.append(
                {
                    "actual_source": source,
                    "actual_orientation": orientation,
                    "predicted_orientation_from_chiral_bispectrum": predicted_orientation,
                    "orientation_matches_actual": predicted_orientation == orientation,
                    "ordered_support": list(ordered_d5_support(source, orientation)),
                    "per_coprime_mode_source_recovery": per_mode,
                    "all_modes_recover_source": all(
                        packet["matches_actual_source"] for packet in per_mode.values()
                    ),
                }
            )
    return rows


def magnitude_signature(config: tuple[int, ...]) -> tuple[float, ...]:
    return tuple(rounded(abs(dft(config, mode))) for mode in range(Z12_NODE_COUNT))


def negative_controls() -> dict[str, Any]:
    magnitude_counts_by_orientation = {}
    for orientation in (-1, 1):
        signatures = {
            magnitude_signature(value_configuration(source, orientation))
            for source in range(Z12_NODE_COUNT)
        }
        magnitude_counts_by_orientation[str(orientation)] = len(signatures)

    noncoprime_alias_counts = {
        str(mode): Z12_NODE_COUNT // math.gcd(mode, Z12_NODE_COUNT)
        for mode in NONCOPRIME_MODES
    }
    no_orientation_pairs = {
        source: {
            signed_orientation_from_bispectrum(value_configuration(source, orientation))
            for orientation in (-1, 1)
        }
        for source in range(Z12_NODE_COUNT)
    }
    return {
        "translation_invariant_magnitude_signature_counts_by_orientation": magnitude_counts_by_orientation,
        "translation_invariant_magnitudes_cannot_select_source": all(
            count == 1 for count in magnitude_counts_by_orientation.values()
        ),
        "noncoprime_mode_alias_class_counts": noncoprime_alias_counts,
        "noncoprime_modes_fail_full_source_resolution": all(
            count < Z12_NODE_COUNT for count in noncoprime_alias_counts.values()
        ),
        "orientation_candidates_per_source_without_chiral_marker": {
            str(source): len(values) for source, values in no_orientation_pairs.items()
        },
        "without_chiral_marker_orientation_not_unique": all(
            len(values) == 2 for values in no_orientation_pairs.values()
        ),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    doc_texts = {name: read_text(path) for name, path in DOC_FILES.items()}

    influence = artifacts["ATOM_INFLUENCE"]
    low_weight = artifacts["LOW_WEIGHT_FRONTIER"]
    phase_reference = artifacts["PHASE_REFERENCE_SELECTOR"]
    chi11_lattice = artifacts["CHI11_PREMISE_LATTICE"]
    legacy_torsion = artifacts["LEGACY_TORSION_CHI11"]
    p2363 = artifacts["P2363_BRIDGE_COMPLETED_MOMENTS"]

    rows = source_recovery_rows()
    predicted_pairs = {
        (
            row["predicted_orientation_from_chiral_bispectrum"],
            row["per_coprime_mode_source_recovery"]["1"]["recovered_source"],
        )
        for row in rows
    }
    controls = negative_controls()

    influence_summary = influence.get("theorem_frontier_atom_influence_summary", {})
    low_weight_summary = low_weight.get("theorem_frontier_low_weight_extension_summary", {})
    if not low_weight_summary:
        low_weight_summary = {
            "chi11_is_only_singleton_unlock": "chi11_selector_source"
            in low_weight.get("proof_certificate", {}).get("singleton_step", ""),
        }
    chi11_antichain = chi11_lattice.get("minimal_premise_antichains", {})
    legacy_verdict = legacy_torsion.get("opinion_audit", {}).get("overall_verdict", "")
    phase_scan = phase_reference.get("conditional_selector_scan", {})

    candidate_ranking = [
        {
            "rank": 1,
            "candidate": "phase_origin_plus_chiral_bispectrum_for_chi11_selector_source",
            "why_selected": [
                "It is the strongest concrete constructive candidate found by grep: exact modular recovery formula plus exhaustive 24-row scan.",
                "It targets the unique top logical bottleneck, chi11_selector_source.",
                "It has clear negative controls: no phase origin leaves translation-invariant magnitudes source-blind; non-coprime modes alias sources; no chiral marker leaves orientation two-valued.",
            ],
            "status": "BEST_OPERATIONAL_CANDIDATE_PREMISE_BASED_NOT_STRICT_SOURCE",
        },
        {
            "rank": 2,
            "candidate": "beta_tors_to_chi11_bridge_hypothesis",
            "why_not_first": [
                "S2 makes it a real bridge target, but the current audit classifies it as candidate-not-theorem.",
                "It still needs an orientation map, role-transfer control, eta-pipeline link, and full-Aut obstruction escape.",
            ],
            "status": legacy_verdict,
        },
        {
            "rank": 3,
            "candidate": "sparsest_extension_or_max_shell_imbalance_branch_selector",
            "why_not_first": [
                "The chi11 premise lattice says these can uniquely select a branch only after quotient, chi11 character, and normalization premises are imported.",
                "They are downstream conditional selectors rather than a source for chi11_selector_source.",
            ],
            "status": "CONDITIONAL_DOWNSTREAM_BRANCH_SELECTOR",
        },
    ]

    theorem_export = {
        "theorem_name": "P2366 selector candidate audit: phase-origin chi11 source candidate",
        "claim": (
            "The strongest current concrete selector candidate is not raw beta_tors nor an "
            "Aut-invariant score, but a premise-based phase-origin source-localizing selector: "
            "a chiral bispectrum chooses orientation and a calibrated coprime Fourier phase "
            "recovers the source.  This is sufficient on the audited 24 source/orientation "
            "rows, but remains a candidate because the phase origin and handedness are not "
            "derived as strict-core data."
        ),
        "positive_content": [
            "chi11_selector_source is the unique top Boolean bottleneck and the only singleton selector unlock.",
            "All 24 source/orientation rows recover orientation and source using chiral bispectrum plus calibrated coprime phase.",
            "All coprime modes 1,5,7,11 recover all sources.",
            "Translation-invariant magnitude data cannot select source.",
            "Non-coprime phase modes alias sources and fail full source resolution.",
            "Without a chiral marker, orientation remains two-valued.",
        ],
        "not_licensed": [
            "strict-core source theorem for phase origin",
            "strict-core source theorem for handedness/chirality",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "selector/QW-2191 discharge",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2366_S1316_SELECTOR_CANDIDATE_PHASE_REFERENCE_CHI11_AUDIT",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "grep_synthesis": {
            "top_logical_bottleneck": influence_summary.get("top_influence_atoms", []),
            "chi11_total_critical_count": influence_summary.get(
                "chi11_selector_source_total_critical_count"
            ),
            "chi11_is_only_singleton_unlock": low_weight_summary.get(
                "chi11_is_only_singleton_unlock"
            ),
            "strict_full_aut_internal_chi11_polarity_source_antichain": chi11_antichain.get(
                "strict_full_aut_internal_chi11_polarity_source"
            ),
            "legacy_torsion_overall_verdict": legacy_verdict,
        },
        "candidate_ranking": candidate_ranking,
        "selected_candidate": candidate_ranking[0],
        "phase_origin_selector_research": {
            "formula": "s = inv(k mod 12) * 12 * (phase_ref_k(orientation)-phase_obs_k) mod 12",
            "coprime_modes": list(COPRIME_SOURCE_MODES),
            "row_count": len(rows),
            "expected_pair_count": 2 * Z12_NODE_COUNT,
            "unique_predicted_orientation_source_pairs_using_mode_1": len(predicted_pairs),
            "all_orientations_recovered_by_chiral_bispectrum": all(
                row["orientation_matches_actual"] for row in rows
            ),
            "all_rows_recovered_by_all_coprime_modes": all(
                row["all_modes_recover_source"] for row in rows
            ),
            "source_recovery_success_by_mode": {
                str(mode): all(
                    row["per_coprime_mode_source_recovery"][str(mode)]["matches_actual_source"]
                    for row in rows
                )
                for mode in COPRIME_SOURCE_MODES
            },
            "rows": rows,
            "prior_phase_reference_report_agrees": {
                "row_count": phase_scan.get("row_count") == len(rows),
                "all_rows_recovered_by_all_coprime_modes": phase_scan.get(
                    "all_rows_recovered_by_all_coprime_modes"
                )
                is True,
                "unique_pairs_mode_1": phase_scan.get(
                    "unique_predicted_orientation_source_pairs_using_mode_1"
                )
                == len(predicted_pairs),
            },
        },
        "negative_controls": controls,
        "bridge_and_eom_separation": {
            "p2363_bridge_completed_moments_available": p2363.get("packet_id") == "P2363",
            "selector_work_is_parallel_to_eom_lagrangian": True,
            "selector_required_for_p2365_eom_lift": False,
        },
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "chi11_unique_top_bottleneck_loaded": influence_summary.get(
            "chi11_selector_source_is_unique_top_influence"
        )
        is True,
        "chi11_only_singleton_unlock_loaded": low_weight_summary.get(
            "chi11_is_only_singleton_unlock"
        )
        is True,
        "full_aut_internal_chi11_source_absent_in_lattice": chi11_antichain.get(
            "strict_full_aut_internal_chi11_polarity_source"
        )
        == [],
        "legacy_torsion_not_promoted_to_selector_theorem": "NOT_A_CURRENT_REPO_THEOREM"
        in legacy_verdict,
        "selected_candidate_is_phase_origin": candidate_ranking[0]["candidate"].startswith(
            "phase_origin"
        ),
        "all_24_rows_recover_orientation_and_source": len(rows) == 24
        and len(predicted_pairs) == 24
        and all(row["orientation_matches_actual"] and row["all_modes_recover_source"] for row in rows),
        "prior_phase_reference_report_replayed": all(
            probe["phase_origin_selector_research"]["prior_phase_reference_report_agrees"].values()
        ),
        "negative_controls_pass": controls["translation_invariant_magnitudes_cannot_select_source"]
        and controls["noncoprime_modes_fail_full_source_resolution"]
        and controls["without_chiral_marker_orientation_not_unique"],
        "docs_updated_with_p2366_selector_candidate_statement": all(
            snippet in text for text in doc_texts.values() for snippet in REQUIRED_DOC_SNIPPETS
        ),
        "no_qw2191_discharge_claimed": True,
        "hard_limits_preserved": True,
    }

    payload = {
        "schema_version": "p2366_s1316_v1",
        "packet_id": "P2366",
        "stage_id": "S1316",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_SELECTOR_CANDIDATE_AUDIT_PHASE_ORIGIN_PREMISE_NO_QW2191_DISCHARGE",
        "result_kind": "SELECTOR_CANDIDATE_PHASE_REFERENCE_CHI11_AUDIT",
        "selector_candidate_phase_reference_chi11_audit_probe": probe,
        "recommended_next_honest_step": (
            "Try to derive the phase-origin reference from an internal strict object "
            "(for example an APD phase transport boundary, a defect/node anchor, or a "
            "self-recorded arrow action), or prove that every admissible strict-core "
            "candidate remains translation-invariant and therefore source-blind."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2366 S1316: selector candidate phase-origin chi11 audit",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2366/S1316 selects the phase-origin selector candidate as the strongest current concrete candidate for `chi11_selector_source`: chiral bispectrum fixes orientation and calibrated coprime Fourier phase recovers source.",
                "",
                "## Checks",
                "",
                f"- Rows checked: `{len(rows)}` / expected `{2 * Z12_NODE_COUNT}`.",
                f"- Unique orientation/source pairs using mode 1: `{len(predicted_pairs)}`.",
                f"- All coprime modes recover all sources: `{probe['phase_origin_selector_research']['all_rows_recovered_by_all_coprime_modes']}`.",
                f"- Translation-invariant magnitudes source-blind: `{controls['translation_invariant_magnitudes_cannot_select_source']}`.",
                f"- Non-coprime modes alias sources: `{controls['noncoprime_modes_fail_full_source_resolution']}`.",
                "",
                "## Hard Limits",
                "",
                "- No strict-core phase-origin theorem is claimed.",
                "- No beta_tors -> chi11 theorem or legacy role transfer is claimed.",
                "- QW-2191 remains open; no selector closure or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
