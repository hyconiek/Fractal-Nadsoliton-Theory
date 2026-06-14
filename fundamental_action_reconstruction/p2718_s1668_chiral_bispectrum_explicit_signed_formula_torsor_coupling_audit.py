#!/usr/bin/env python3
"""P2718/S1668: chiral-bispectrum explicit signed formula torsor-coupling audit.

P2717 said not to keep listing generic pseudoscalar names; the next admissible
move needs one explicit signed formula/artifact.  P2718 uses the existing
phase-origin/chiral-bispectrum finite formula family as that concrete formula,
recomputes its signed marker on the 24 source/orientation rows, and then checks
whether the marker is already a strict non-premise source coupled to the
P2708/P2714 orientation torsor.
"""
from __future__ import annotations

import cmath
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.json"
MD = GEN / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
DISTANCE_SELECTED = 5
FORWARD_ASSIGNMENT = (2, 2, 2, 1, 1)
ORIENTATION_BISPECTRUM_PAIR = (1, 5)
ROUND_DIGITS = 9

INPUTS = {
    "P2717_SOURCE_MATRIX": GEN / "p2717_s1667_concrete_pseudoscalar_chiral_source_candidate_matrix.json",
    "P2716_PSEUDOSCALAR_ACCEPTANCE": GEN / "p2716_s1666_inversion_odd_pseudoscalar_source_acceptance_audit.json",
    "P2412_CHI11_SCOPE": GEN / "p2412_s1362_chi11_selector_scope_separation_certificate.json",
    "P2367_PHASE_ORIGIN_NO_GO": GEN / "p2367_s1317_selector_phase_origin_admissibility_no_go_probe.json",
    "P2708_BOUNDARY_COCYCLE": GEN / "p2708_s1658_z12_boundary_1_cocycle_selector_source_obstruction.json",
    "P2714_ORIENTATION_TORSOR": GEN / "p2714_s1664_z12_orientation_torsor_global_section_obstruction.json",
}

ACCEPTANCE_CRITERIA = [
    "explicit_signed_formula_computed",
    "nonzero_signed_value_on_all_rows",
    "translation_or_source_localizer_strictly_exported",
    "phase_origin_reference_nonpremise",
    "coupling_to_p2708_p2714_torsor_exported",
    "qw2191_safe_nonpremise_selector_source",
]

NEGATIVE_EXPORT_FLAGS = [
    "chiral_bispectrum_accepted_as_strict_source",
    "translation_or_source_localizer_strictly_exported",
    "phase_origin_reference_nonpremise",
    "coupling_to_p2708_p2714_torsor_exported",
    "strict_mechanism_fixing_lambda_exported",
    "qw2191_discharged",
    "pair12_strict_core_upgrade_exported",
    "bridge_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


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
        config, (left_mode + right_mode) % Z12_NODE_COUNT
    ).conjugate()


def marker_rows() -> list[dict[str, Any]]:
    rows = []
    for orientation in (-1, 1):
        for source in range(Z12_NODE_COUNT):
            marker = bispectrum(value_configuration(source, orientation), *ORIENTATION_BISPECTRUM_PAIR)
            signed_value = rounded(marker.imag)
            rows.append({
                "source": source,
                "orientation": orientation,
                "bispectrum_pair": list(ORIENTATION_BISPECTRUM_PAIR),
                "marker_real": rounded(marker.real),
                "marker_imag": signed_value,
                "nonzero_signed_marker": abs(signed_value) > 1e-9,
                "orientation_recovered_from_marker": -1 if signed_value > 0 else 1,
            })
    return rows


def translation_degeneracy_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    by_orientation: dict[int, set[float]] = {-1: set(), 1: set()}
    by_source: dict[int, set[float]] = {source: set() for source in range(Z12_NODE_COUNT)}
    for row in rows:
        by_orientation[row["orientation"]].add(row["marker_imag"])
        by_source[row["source"]].add(row["marker_imag"])
    return {
        "orientation_marker_values": {str(k): sorted(v) for k, v in by_orientation.items()},
        "source_marker_value_counts": {str(k): len(v) for k, v in by_source.items()},
        "orientation_separating": by_orientation[-1] == {2.0} and by_orientation[1] == {-2.0},
        "source_localizing": all(len(v) == Z12_NODE_COUNT for v in by_source.values()),
        "source_blind_under_translation": all(len(v) == 2 for v in by_source.values()),
    }


def acceptance_matrix(rows: list[dict[str, Any]], degeneracy: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "explicit_signed_formula_computed": True,
        "nonzero_signed_value_on_all_rows": all(row["nonzero_signed_marker"] for row in rows),
        "translation_or_source_localizer_strictly_exported": False,
        "phase_origin_reference_nonpremise": False,
        "coupling_to_p2708_p2714_torsor_exported": False,
        "qw2191_safe_nonpremise_selector_source": False,
    }
    missing = [criterion for criterion in ACCEPTANCE_CRITERIA if not facts[criterion]]
    return {
        "facts": facts,
        "missing_criteria": missing,
        "accepted_as_strict_pseudoscalar_source": not missing,
        "blocker": "The chiral-bispectrum imaginary marker is a real nonzero signed formula and separates orientation, but it is translation/source-blind without a non-premise phase-origin/source localizer and has no exported coupling theorem to the P2708/P2714 torsor.",
        "degeneracy_used": degeneracy,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2718/S1668 chiral-bispectrum explicit signed formula torsor-coupling audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Formula result",
        f"- checked_rows={payload['finite_summary']['checked_rows']}",
        f"- orientation_separating={payload['finite_summary']['orientation_separating']}",
        f"- nonzero_signed_value_on_all_rows={payload['acceptance_matrix']['facts']['nonzero_signed_value_on_all_rows']}",
        "",
        "## Acceptance",
        f"- accepted_as_strict_pseudoscalar_source={payload['acceptance_matrix']['accepted_as_strict_pseudoscalar_source']}",
        f"- missing={payload['acceptance_matrix']['missing_criteria']}",
        payload['acceptance_matrix']['blocker'],
        "",
        "## Decision",
        payload["decision"]["reason"],
        "",
        "## Next honest step",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    rows = marker_rows()
    degeneracy = translation_degeneracy_summary(rows)
    acceptance = acceptance_matrix(rows, degeneracy)
    no_unlock = (
        len(rows) == 24
        and degeneracy["orientation_separating"]
        and acceptance["facts"]["nonzero_signed_value_on_all_rows"]
        and not acceptance["accepted_as_strict_pseudoscalar_source"]
    )
    payload = {
        "status": "P2718_CHIRAL_BISPECTRUM_SIGNED_FORMULA_POSITIVE_BUT_NO_STRICT_TORSOR_SOURCE" if no_unlock else "P2718_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: read_json(path).get("status") for key, path in INPUTS.items()},
        "finite_summary": {
            "formula": "Im(B_{1,5}) for the P2367/P2366 chiral-bispectrum marker on 12-node D5 supports",
            "checked_rows": len(rows),
            "orientation_separating": degeneracy["orientation_separating"],
            "orientation_marker_values": degeneracy["orientation_marker_values"],
            "source_blind_under_translation": degeneracy["source_blind_under_translation"],
        },
        "marker_rows": rows,
        "translation_degeneracy_summary": degeneracy,
        "acceptance_matrix": acceptance,
        "decision": {
            "explicit_signed_formula_computed": True,
            "chiral_bispectrum_accepted_as_strict_source": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2718 tests the concrete chiral-bispectrum formula demanded after P2717.  Its imaginary marker is nonzero on all 24 source/orientation rows and separates orientation with values +2 and -2.  However, the marker is still not a strict torsor-breaking source: it lacks a non-premise phase-origin/source localizer and no theorem couples it to the P2708/P2714 orientation torsor as a QW-2191-safe selector source.",
            "next_honest_step": "Do not promote the chiral-bispectrum marker by itself.  The next admissible move is a narrow phase-origin/source-localizer theorem audit for this exact formula, proving a non-premise origin reference and exported torsor coupling; if that cannot be supplied, preserve the P2697-P2718 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2718/S1668 chiral-bispectrum explicit signed formula torsor-coupling audit", "## P2718/S1668 chiral-bispectrum explicit signed formula torsor-coupling audit\n\n`P2718/S1668` tests one explicit signed formula after P2717: the P2366/P2367 chiral-bispectrum marker `Im(B_{1,5})` on the 12-node D5 support family.  The marker is nonzero on all 24 source/orientation rows and separates orientation by the values `+2` and `-2`.  This is real signed formula evidence, but it is still not a strict torsor-breaking source because a non-premise phase-origin/source localizer and an exported coupling theorem to the P2708/P2714 orientation torsor remain missing.  No strict `lambda` fixing, `QW-2191` discharge, pair12 strict-core upgrade, bridge closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2718/S1668 chiral-bispectrum formula Ltotal guard", "## P2718/S1668 chiral-bispectrum formula Ltotal guard\n\n`P2718/S1668` recomputes a concrete signed chiral-bispectrum marker, but it is not a variational source construction.  Without a non-premise phase-origin/source localizer and torsor-coupling theorem, the marker does not promote `L_total`, selector closure, pair12 strict-core, role transfer, bridge closure, or ToE.\n")
    append_once(AGENTS, "Current chiral-bispectrum signed formula torsor-coupling guardrail (P2718/S1668, 2026-06-14)", "## Current chiral-bispectrum signed formula torsor-coupling guardrail (P2718/S1668, 2026-06-14)\n\n- P2718 recomputes the explicit P2366/P2367 chiral-bispectrum marker `Im(B_{1,5})`: it is nonzero on all 24 source/orientation rows and separates orientation by `+2/-2`.\n- The marker is not yet a strict torsor-breaking source because it lacks a non-premise phase-origin/source localizer and an exported coupling theorem to the P2708/P2714 orientation torsor; therefore it does not fix `lambda` or discharge `QW-2191`.\n- Do not promote the chiral-bispectrum marker alone to selector closure, pair12 strict-core, role transfer, bridge closure, `L_total`, or ToE; a next admissible move is only a narrow phase-origin/source-localizer theorem audit for this exact formula, or preservation of the P2697-P2718 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
