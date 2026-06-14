#!/usr/bin/env python3
"""P2727/S1677: orientation-transition law equivariance and polarity no-go.

P2726 found a real nonzero +/-4 marker jump, but only when an orientation flip
is imported.  P2727 asks the next proof-grade finite question: can a simple
source-independent orientation-transition law on the two orientation states
select a nonzero polarity without becoming a premise selector?

We enumerate the four possible maps {-,+}->{-,+} (preserve, flip, collapse to
+, collapse to -), combine each with all 12 Z12 source translations, and score
the P2718 marker jump.  The Aut/inversion-equivariant laws are exactly preserve
and flip: preserve gives zero; flip gives balanced +/-4 and selects no unique
polarity.  The two polarity-selecting collapse laws fail inversion equivariance,
so they are explicit premises rather than strict sources.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any, Callable

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2727_s1677_orientation_transition_law_equivariance_and_polarity_no_go.json"
MD = GEN / "p2727_s1677_orientation_transition_law_equivariance_and_polarity_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

Z12 = 12
ORIENTATIONS = (-1, 1)
VELOCITIES = tuple(range(Z12))
INPUTS = {
    "P2726_AFFINE_ORIENTATION_FLOW": GEN / "p2726_s1676_chiral_bispectrum_affine_orientation_flow_transition_matrix.json",
    "P2721_SIGN_TORSOR_COUPLING": GEN / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.json",
    "P2718_CHIRAL_BISPECTRUM_MARKER": GEN / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "equivariant_polarity_selecting_orientation_law_exported",
    "nonpremise_orientation_branch_selected",
    "p2721_coupling_polarity_selected",
    "strict_mechanism_fixing_lambda_exported",
    "qw2191_discharged",
    "pair12_strict_core_upgrade_exported",
    "bridge_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]

LAW_TABLE: dict[str, dict[int, int]] = {
    "preserve_orientation": {-1: -1, 1: 1},
    "flip_orientation": {-1: 1, 1: -1},
    "collapse_to_plus_orientation": {-1: 1, 1: 1},
    "collapse_to_minus_orientation": {-1: -1, 1: -1},
}


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"missing": True, "path": rel(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def marker_table(p2718: dict[str, Any]) -> dict[tuple[int, int], float]:
    return {
        (int(row["orientation"]), int(row["source"])): float(row["marker_imag"])
        for row in p2718.get("marker_rows", [])
    }


def inversion_equivariant(mapping: dict[int, int]) -> bool:
    # Inversion sends orientation o -> -o.  Equivariance requires f(-o) = -f(o).
    return all(mapping[-orientation] == -mapping[orientation] for orientation in ORIENTATIONS)


def law_rows(table: dict[tuple[int, int], float]) -> dict[str, Any]:
    rows = []
    law_summaries = []
    for law_name, mapping in LAW_TABLE.items():
        deltas = []
        for velocity in VELOCITIES:
            for orientation in ORIENTATIONS:
                for source in range(Z12):
                    next_orientation = mapping[orientation]
                    next_source = (source + velocity) % Z12
                    delta = round(table[(next_orientation, next_source)] - table[(orientation, source)], 9)
                    if delta == -0.0:
                        delta = 0.0
                    deltas.append(delta)
                    rows.append({
                        "law": law_name,
                        "velocity": velocity,
                        "orientation": orientation,
                        "source": source,
                        "next_orientation": next_orientation,
                        "next_source": next_source,
                        "delta_marker_imag": delta,
                    })
        histogram = dict(sorted(Counter(deltas).items()))
        nonzero = [delta for delta in deltas if abs(delta) > 1e-9]
        positive = [delta for delta in deltas if delta > 0]
        negative = [delta for delta in deltas if delta < 0]
        law_summaries.append({
            "law": law_name,
            "mapping": {str(k): v for k, v in mapping.items()},
            "inversion_equivariant": inversion_equivariant(mapping),
            "row_count": len(deltas),
            "delta_values": sorted(set(deltas)),
            "delta_histogram": histogram,
            "nonzero_delta_count": len(nonzero),
            "positive_delta_count": len(positive),
            "negative_delta_count": len(negative),
            "balanced_nonzero_polarities": len(positive) == len(negative) and len(positive) > 0,
            "selects_single_nonzero_polarity": (len(positive) > 0 and len(negative) == 0) or (len(negative) > 0 and len(positive) == 0),
        })
    equivariant = [row for row in law_summaries if row["inversion_equivariant"]]
    equivariant_selectors = [row for row in equivariant if row["selects_single_nonzero_polarity"]]
    premise_selectors = [row for row in law_summaries if row["selects_single_nonzero_polarity"] and not row["inversion_equivariant"]]
    return {
        "checked_transition_rows": len(rows),
        "law_count": len(LAW_TABLE),
        "velocity_count": len(VELOCITIES),
        "rows_per_law": len(VELOCITIES) * len(ORIENTATIONS) * Z12,
        "law_summaries": law_summaries,
        "equivariant_law_names": [row["law"] for row in equivariant],
        "equivariant_polarity_selecting_law_count": len(equivariant_selectors),
        "premise_polarity_selecting_law_names": [row["law"] for row in premise_selectors],
        "obstruction": "The only inversion-equivariant source-independent laws are preserve and flip.  Preserve has zero jump; flip has balanced +4/-4 jumps.  The laws that select one nonzero polarity collapse to a chosen orientation and fail inversion equivariance, so they are premise selectors.",
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "all_source_independent_orientation_laws_enumerated": audit["law_count"] == 4,
        "all_z12_velocities_checked": audit["velocity_count"] == 12,
        "equivariant_laws_identified": audit["equivariant_law_names"] == ["preserve_orientation", "flip_orientation"],
        "equivariant_polarity_selecting_law_exported": audit["equivariant_polarity_selecting_law_count"] > 0,
        "nonpremise_orientation_branch_selected": False,
        "p2721_coupling_polarity_selected": False,
    }
    missing = [criterion for criterion, value in facts.items() if not value]
    return {
        "facts": facts,
        "missing_criteria": missing,
        "accepted_as_strict_orientation_law": not missing,
        "blocker": "Finite enumeration finds no inversion-equivariant source-independent orientation-transition law that selects a single nonzero marker-jump polarity.  Polarity-selecting laws are exactly orientation-collapse premises.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["orientation_law_audit"]
    lines = [
        "# P2727/S1677 orientation-transition law equivariance and polarity no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite law enumeration",
        f"- checked_transition_rows={audit['checked_transition_rows']}",
        f"- law_count={audit['law_count']}",
        f"- equivariant_law_names={audit['equivariant_law_names']}",
        f"- equivariant_polarity_selecting_law_count={audit['equivariant_polarity_selecting_law_count']}",
        f"- premise_polarity_selecting_law_names={audit['premise_polarity_selecting_law_names']}",
        audit["obstruction"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    audit = law_rows(marker_table(inputs["P2718_CHIRAL_BISPECTRUM_MARKER"]))
    acceptance = acceptance_matrix(audit)
    no_go = audit["equivariant_polarity_selecting_law_count"] == 0 and not acceptance["accepted_as_strict_orientation_law"]
    payload = {
        "status": "P2727_ORIENTATION_TRANSITION_LAW_EQUIVARIANCE_POLARITY_NO_GO" if no_go else "P2727_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "all source-independent orientation-transition laws on {-1,+1}, crossed with all Z12 source velocities",
        "orientation_law_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "equivariant_polarity_selecting_orientation_law_exported": False,
            "nonpremise_orientation_branch_selected": False,
            "p2721_coupling_polarity_selected": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2727 exhausts the finite source-independent orientation-law class.  The only inversion-equivariant laws are preserve and flip: preserve has zero marker jump, while flip has balanced +4/-4 jumps and no single polarity.  The only single-polarity laws collapse to a chosen orientation, hence import the missing selector premise.",
            "next_honest_step": "Do not continue source-independent orientation-law variants.  A next admissible proof-grade move must introduce a source-dependent but non-premise strict invariant that breaks inversion equivariance with a computable signed value and an exported P2721 polarity-coupling theorem; otherwise preserve the P2697-P2727 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2727/S1677 orientation-transition law equivariance and polarity no-go", "## P2727/S1677 orientation-transition law equivariance and polarity no-go\n\n`P2727/S1677` exhausts source-independent orientation-transition laws on `{-1,+1}` crossed with all Z12 source velocities.  The inversion-equivariant laws are only preserve and flip: preserve gives zero marker jump, while flip gives balanced `+4/-4` jumps and no unique polarity.  The single-polarity laws collapse to a chosen orientation and fail inversion equivariance, so they are premise selectors.  No P2721 polarity selection, `lambda` fixing, `QW-2191` discharge, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2727/S1677 no source-independent orientation law Ltotal guard", "## P2727/S1677 no source-independent orientation law Ltotal guard\n\n`P2727/S1677` finds no inversion-equivariant source-independent orientation-transition law that selects one nonzero marker-jump polarity.  Since the polarity-selecting laws are premise selectors, they cannot be promoted to a variational source term or `L_total`.\n")
    append_once(AGENTS, "Current orientation-transition law equivariance no-go guardrail (P2727/S1677, 2026-06-14)", "## Current orientation-transition law equivariance no-go guardrail (P2727/S1677, 2026-06-14)\n\n- P2727 exhausts source-independent orientation-transition laws on `{-1,+1}` crossed with all Z12 source velocities (1152 rows total).\n- The only inversion-equivariant laws are preserve and flip: preserve has zero marker jump, and flip has balanced `+4/-4` jumps; the laws selecting one nonzero polarity collapse to a chosen orientation and fail inversion equivariance.\n- Do not continue source-independent orientation-law variants as `QW-2191` discharge, selector closure, pair12 strict-core upgrade, role transfer, bridge closure, `L_total`, or ToE.  A next admissible move must introduce a source-dependent but non-premise strict invariant with a computable signed value and an exported P2721 polarity-coupling theorem, or preserve the P2697-P2727 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
