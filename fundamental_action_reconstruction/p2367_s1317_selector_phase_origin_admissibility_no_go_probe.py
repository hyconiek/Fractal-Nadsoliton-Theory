#!/usr/bin/env python3
from __future__ import annotations

import cmath
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
SCRATCH = ROOT / "scratch"

OUT = GEN / "p2367_s1317_selector_phase_origin_admissibility_no_go_probe.json"
MD = GEN / "p2367_s1317_selector_phase_origin_admissibility_no_go_probe.md"

SOURCE_FILES = {
    "P2366_SELECTOR_CANDIDATE": GEN
    / "p2366_s1316_selector_candidate_phase_reference_chi11_audit_probe.json",
    "C12_SOURCE_TRANSLATION_NO_GO": SCRATCH
    / "bridge_strict_alpha_c12_source_translation_no_go_report.json",
    "PHASE_REFERENCE_SELECTOR": SCRATCH
    / "bridge_strict_alpha_phase_reference_source_selector_certificate_report.json",
    "CHI11_PREMISE_LATTICE": SCRATCH
    / "bridge_strict_alpha_chi11_premise_dependency_lattice_audit_report.json",
}

DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
DISTANCE_SELECTED = 5
FORWARD_ASSIGNMENT = (2, 2, 2, 1, 1)
COPRIME_SOURCE_MODES = (1, 5, 7, 11)
NONCOPRIME_MODES = (2, 3, 4, 6, 8, 9, 10)
ORIENTATION_BISPECTRUM_PAIR = (1, 5)
ROUND_DIGITS = 9


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode(
            "utf-8"
        )
    ).hexdigest()


def rounded(value: float) -> float:
    result = round(value, ROUND_DIGITS)
    return 0.0 if result == -0.0 else result


def rounded_complex(value: complex) -> tuple[float, float]:
    return (rounded(value.real), rounded(value.imag))


def ordered_d5_support(source: int, orientation: int) -> tuple[int, ...]:
    step = (orientation * DISTANCE_SELECTED) % Z12_NODE_COUNT
    return tuple((source + index * step) % Z12_NODE_COUNT for index in range(SUPPORT_SIZE))


def value_configuration(source: int, orientation: int) -> tuple[int, ...]:
    values = [0] * Z12_NODE_COUNT
    for node, value in zip(ordered_d5_support(source, orientation), FORWARD_ASSIGNMENT):
        values[node] = value
    return tuple(values)


def translate_config(config: tuple[int, ...], shift: int) -> tuple[int, ...]:
    translated = [0] * Z12_NODE_COUNT
    for node, value in enumerate(config):
        translated[(node + shift) % Z12_NODE_COUNT] = value
    return tuple(translated)


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


def power_signature(config: tuple[int, ...]) -> tuple[float, ...]:
    return tuple(rounded(abs(dft(config, mode)) ** 2) for mode in range(Z12_NODE_COUNT))


def complete_bispectrum_signature(config: tuple[int, ...]) -> tuple[tuple[float, float], ...]:
    return tuple(
        rounded_complex(bispectrum(config, left_mode, right_mode))
        for left_mode in range(Z12_NODE_COUNT)
        for right_mode in range(Z12_NODE_COUNT)
    )


def phase_signature(config: tuple[int, ...], modes: Iterable[int]) -> tuple[float, ...]:
    return tuple(phase_turns(dft(config, mode)) for mode in modes)


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
    return int(round(((reference_phase - phase) % 1.0) * Z12_NODE_COUNT * inverse)) % Z12_NODE_COUNT


def orbit_certificate() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for orientation in (-1, 1):
        base = value_configuration(0, orientation)
        translated_sources = []
        stabilizer = []
        for shift in range(Z12_NODE_COUNT):
            translated = translate_config(base, shift)
            if translated == base:
                stabilizer.append(shift)
            for source in range(Z12_NODE_COUNT):
                if translated == value_configuration(source, orientation):
                    translated_sources.append(source)
        rows.append(
            {
                "orientation": orientation,
                "orbit_size_from_source_rows": len(
                    {value_configuration(source, orientation) for source in range(Z12_NODE_COUNT)}
                ),
                "translated_sources_from_source_0": translated_sources,
                "translation_action_is_transitive_on_sources": sorted(translated_sources)
                == list(range(Z12_NODE_COUNT)),
                "stabilizer_of_source_0": stabilizer,
                "translation_action_is_free_on_sources": stabilizer == [0],
            }
        )
    return rows


def admissibility_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    references = calibrated_phase_references()
    for orientation in (-1, 1):
        configs = [value_configuration(source, orientation) for source in range(Z12_NODE_COUNT)]
        powers = {power_signature(config) for config in configs}
        bispectra = {complete_bispectrum_signature(config) for config in configs}
        chiral_signs = {signed_orientation_from_bispectrum(config) for config in configs}
        coprime_phase_signatures = {phase_signature(config, COPRIME_SOURCE_MODES) for config in configs}
        noncoprime_alias_counts = {
            str(mode): len({phase_signature(config, (mode,)) for config in configs})
            for mode in NONCOPRIME_MODES
        }
        recovered_by_mode = defaultdict(list)
        for source, config in enumerate(configs):
            predicted_orientation = signed_orientation_from_bispectrum(config)
            for mode in COPRIME_SOURCE_MODES:
                observed = phase_turns(dft(config, mode))
                reference = references[str(predicted_orientation)][str(mode)]
                recovered_by_mode[str(mode)].append(recover_source_from_phase(observed, reference, mode))
        rows.append(
            {
                "orientation": orientation,
                "unique_power_signatures": len(powers),
                "unique_complete_bispectrum_signatures": len(bispectra),
                "unique_chiral_bispectrum_orientations": sorted(chiral_signs),
                "unique_raw_coprime_phase_signatures": len(coprime_phase_signatures),
                "noncoprime_phase_alias_class_counts": noncoprime_alias_counts,
                "calibrated_coprime_phase_recovers_all_sources_by_mode": {
                    mode: sorted(sources) == list(range(Z12_NODE_COUNT))
                    for mode, sources in recovered_by_mode.items()
                },
            }
        )
    return rows


def max_covariance_errors() -> dict[str, float]:
    max_power_error = 0.0
    max_bispectrum_error = 0.0
    max_dft_covariance_error = 0.0
    for orientation in (-1, 1):
        for source in range(Z12_NODE_COUNT):
            config = value_configuration(source, orientation)
            for shift in range(Z12_NODE_COUNT):
                translated = translate_config(config, shift)
                for mode in range(Z12_NODE_COUNT):
                    power_error = abs(abs(dft(translated, mode)) ** 2 - abs(dft(config, mode)) ** 2)
                    max_power_error = max(max_power_error, power_error)
                    expected_dft = cmath.exp(-2j * math.pi * mode * shift / Z12_NODE_COUNT) * dft(
                        config,
                        mode,
                    )
                    max_dft_covariance_error = max(
                        max_dft_covariance_error,
                        abs(dft(translated, mode) - expected_dft),
                    )
                    for right_mode in range(Z12_NODE_COUNT):
                        max_bispectrum_error = max(
                            max_bispectrum_error,
                            abs(
                                bispectrum(translated, mode, right_mode)
                                - bispectrum(config, mode, right_mode)
                            ),
                        )
    return {
        "max_power_invariance_error": rounded(max_power_error),
        "max_complete_bispectrum_invariance_error": rounded(max_bispectrum_error),
        "max_dft_translation_covariance_error": rounded(max_dft_covariance_error),
    }


def main() -> None:
    GEN.mkdir(exist_ok=True)
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    docs = {name: path.read_text(encoding="utf-8") if path.exists() else "" for name, path in DOC_FILES.items()}

    orbit_rows = orbit_certificate()
    admissibility = admissibility_rows()
    errors = max_covariance_errors()

    theorem_export = {
        "theorem_name": "P2367 finite source-localization admissibility no-go boundary",
        "claim": (
            "On the audited Z12 selector support, C12 translations act freely and transitively "
            "on the twelve source choices for each fixed chiral orientation. Therefore any "
            "translation-invariant strict-core score is source-blind. Complete Fourier powers "
            "and complete bispectra are verified examples of this no-go; raw coprime phases "
            "can localize source only after an external or internally derived phase-origin "
            "reference is supplied."
        ),
        "positive_content": [
            "The source orbit is a free transitive C12 orbit for both orientations.",
            "Complete power spectra and complete bispectrum signatures are invariant and source-blind.",
            "The chiral bispectrum keeps the orientation information but does not supply a source origin.",
            "Non-coprime phase modes alias source classes.",
            "Coprime phases recover all sources only with a calibrated phase-origin reference.",
        ],
        "not_licensed": [
            "strict-core phase-origin theorem",
            "strict-core selector closure",
            "QW-2191 discharge",
            "beta_tors -> chi11 theorem",
            "legacy physical-role transfer",
            "ToE closure",
        ],
    }

    probe = {
        "probe_id": "P2367_S1317_SELECTOR_PHASE_ORIGIN_ADMISSIBILITY_NO_GO",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "docs": {name: rel(path) for name, path in DOC_FILES.items()},
        "orbit_certificate": orbit_rows,
        "admissibility_rows": admissibility,
        "covariance_and_invariance_errors": errors,
        "candidate_boundary_table": [
            {
                "candidate_class": "translation_invariant_power_or_histogram_score",
                "verdict": "REJECT_AS_SOURCE_LOCALIZER",
                "reason": "constant on each free transitive source orbit",
            },
            {
                "candidate_class": "complete_translation_invariant_bispectrum_score",
                "verdict": "REJECT_AS_SOURCE_LOCALIZER_BUT_KEEP_AS_CHIRAL_ORIENTATION_MARKER",
                "reason": "bispectrum is translation-invariant and source-blind, although its chiral sign fixes orientation",
            },
            {
                "candidate_class": "noncoprime_raw_phase_mode",
                "verdict": "REJECT_AS_FULL_SOURCE_LOCALIZER",
                "reason": "non-coprime modes have nontrivial gcd aliases on Z12",
            },
            {
                "candidate_class": "coprime_phase_with_phase_origin_reference",
                "verdict": "BEST_OPERATIONAL_CANDIDATE_BUT_PREMISE_BASED",
                "reason": "recovers all sources, but imports the still-unproven phase-origin reference",
            },
        ],
        "theorem_export": theorem_export,
    }
    probe["theorem_fingerprint_sha256"] = sha256_json(theorem_export)

    gatekeeper_checks = {
        "source_artifacts_present": all("_missing" not in artifact for artifact in artifacts.values()),
        "p2366_was_run_and_loaded": artifacts["P2366_SELECTOR_CANDIDATE"].get("packet_id") == "P2366",
        "free_transitive_source_orbits": all(
            row["translation_action_is_transitive_on_sources"]
            and row["translation_action_is_free_on_sources"]
            for row in orbit_rows
        ),
        "translation_invariant_examples_source_blind": all(
            row["unique_power_signatures"] == 1 and row["unique_complete_bispectrum_signatures"] == 1
            for row in admissibility
        ),
        "chiral_marker_keeps_orientation": all(
            row["unique_chiral_bispectrum_orientations"] == [row["orientation"]]
            for row in admissibility
        ),
        "noncoprime_modes_alias": all(
            all(count < Z12_NODE_COUNT for count in row["noncoprime_phase_alias_class_counts"].values())
            for row in admissibility
        ),
        "coprime_phase_reference_recovers_sources": all(
            all(row["calibrated_coprime_phase_recovers_all_sources_by_mode"].values())
            for row in admissibility
        ),
        "covariance_errors_zero": all(value == 0.0 for value in errors.values()),
        "docs_updated_with_p2367_boundary": all("P2367/S1317" in text for text in docs.values()),
        "no_qw2191_discharge_claimed": True,
    }

    payload = {
        "schema_version": "p2367_s1317_v1",
        "packet_id": "P2367",
        "stage_id": "S1317",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-06-01T00:00:00+00:00",
        "status": "OPEN_PROGRESS_SELECTOR_ADMISSIBILITY_NO_GO_BOUNDARY_NO_QW2191_DISCHARGE",
        "result_kind": "SELECTOR_PHASE_ORIGIN_ADMISSIBILITY_NO_GO_BOUNDARY",
        "selector_phase_origin_admissibility_no_go_probe": probe,
        "recommended_next_honest_step": (
            "Search for, or refute within a named admissible class, an internal strict phase-origin "
            "provider: APD phase-transport boundary, defect/node anchor, self-recorded arrow action, "
            "or another non-translation-invariant object that is not a silent legacy-role transfer."
        ),
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_PROGRESS_WITH_TRACE_NO_SELECTOR_OR_TOE_CLOSURE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2367 S1317: selector phase-origin admissibility no-go boundary",
                "",
                f"Status: `{payload['status']}`",
                "",
                "## Result",
                "",
                "P2367/S1317 continues P2366 by proving the finite source-localization boundary: for each fixed chiral orientation, the twelve source choices form one free transitive `C12` translation orbit.  Any translation-invariant strict-core score is therefore source-blind.",
                "",
                "## Computed boundary",
                "",
                f"- Free/transitive source orbits: `{gatekeeper_checks['free_transitive_source_orbits']}`.",
                f"- Complete power/bispectrum invariants source-blind: `{gatekeeper_checks['translation_invariant_examples_source_blind']}`.",
                f"- Chiral bispectrum keeps orientation: `{gatekeeper_checks['chiral_marker_keeps_orientation']}`.",
                f"- Non-coprime phase modes alias: `{gatekeeper_checks['noncoprime_modes_alias']}`.",
                f"- Coprime phase with calibrated origin recovers sources: `{gatekeeper_checks['coprime_phase_reference_recovers_sources']}`.",
                "",
                "## Hard limits",
                "",
                "- The result is a no-go boundary plus a positive premise boundary, not a strict-core phase-origin theorem.",
                "- No `beta_tors -> chi11` theorem, legacy role transfer, QW-2191 discharge, selector closure, or ToE closure is claimed.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
