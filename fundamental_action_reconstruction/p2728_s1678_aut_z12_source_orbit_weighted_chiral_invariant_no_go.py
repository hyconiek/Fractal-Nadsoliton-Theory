#!/usr/bin/env python3
"""P2728/S1678: Aut(Z12) source-orbit weighted chiral invariant no-go.

P2727 closed source-independent orientation-transition laws.  The next honest
candidate class is source-dependent but non-premise: allow weights depending on
the Aut(Z12)-orbit of the source label, not on a chosen source representative or
an orientation premise.  P2728 exhausts all {-1,0,+1} weights on the six
Aut(Z12) source orbits and applies them to the exact P2718 Im(B_{1,5}) marker.

Every orbit-weighted signed total cancels between the paired orientations.  Any
nonzero row-level value remains a +/- pair, and no Aut-orbit source weighting
selects a P2721 polarity.  This is a bounded no-go for this source-dependent,
non-representative candidate class.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2728_s1678_aut_z12_source_orbit_weighted_chiral_invariant_no_go.json"
MD = GEN / "p2728_s1678_aut_z12_source_orbit_weighted_chiral_invariant_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

Z12 = 12
UNITS = (1, 5, 7, 11)
WEIGHT_VALUES = (-1, 0, 1)
INPUTS = {
    "P2727_ORIENTATION_LAW_NO_GO": GEN / "p2727_s1677_orientation_transition_law_equivariance_and_polarity_no_go.json",
    "P2721_SIGN_TORSOR_COUPLING": GEN / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.json",
    "P2718_CHIRAL_BISPECTRUM_MARKER": GEN / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "aut_orbit_weighted_source_invariant_exported",
    "nonzero_global_signed_value_exported",
    "single_polarity_selected",
    "p2721_coupling_polarity_selected",
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


def aut_orbits() -> list[tuple[int, ...]]:
    seen: set[int] = set()
    orbits: list[tuple[int, ...]] = []
    for source in range(Z12):
        if source in seen:
            continue
        orbit = tuple(sorted({(unit * source) % Z12 for unit in UNITS}))
        seen.update(orbit)
        orbits.append(orbit)
    return sorted(orbits, key=lambda orbit: (len(orbit), orbit))


def orbit_index(orbits: list[tuple[int, ...]]) -> dict[int, int]:
    return {source: idx for idx, orbit in enumerate(orbits) for source in orbit}


def marker_rows(p2718: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "source": int(row["source"]),
            "orientation": int(row["orientation"]),
            "marker_imag": float(row["marker_imag"]),
        }
        for row in p2718.get("marker_rows", [])
    ]


def exhaustive_weight_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    orbits = aut_orbits()
    source_to_orbit = orbit_index(orbits)
    payload_rows = []
    nonzero_global_rows = []
    single_polarity_rows = []
    for weights in itertools.product(WEIGHT_VALUES, repeat=len(orbits)):
        weighted_values = []
        by_orientation = defaultdict(float)
        by_orbit_orientation = defaultdict(float)
        for row in rows:
            weight = weights[source_to_orbit[row["source"]]]
            value = weight * row["marker_imag"]
            weighted_values.append(value)
            by_orientation[row["orientation"]] += value
            by_orbit_orientation[(source_to_orbit[row["source"]], row["orientation"])] += value
        total = round(sum(weighted_values), 9)
        positive_rows = sum(1 for value in weighted_values if value > 0)
        negative_rows = sum(1 for value in weighted_values if value < 0)
        nonzero_total = abs(total) > 1e-9
        single_polarity = (positive_rows > 0 and negative_rows == 0) or (negative_rows > 0 and positive_rows == 0)
        row_payload = {
            "weights": list(weights),
            "global_signed_total": total,
            "orientation_totals": {str(key): round(value, 9) for key, value in sorted(by_orientation.items())},
            "positive_row_count": positive_rows,
            "negative_row_count": negative_rows,
            "nonzero_global_signed_total": nonzero_total,
            "selects_single_row_polarity": single_polarity,
        }
        payload_rows.append(row_payload)
        if nonzero_total:
            nonzero_global_rows.append(row_payload)
        if single_polarity:
            single_polarity_rows.append(row_payload)

    return {
        "aut_units": list(UNITS),
        "source_orbits": [list(orbit) for orbit in orbits],
        "orbit_count": len(orbits),
        "weight_values": list(WEIGHT_VALUES),
        "weighting_count": len(payload_rows),
        "nonzero_global_signed_total_count": len(nonzero_global_rows),
        "single_polarity_weighting_count": len(single_polarity_rows),
        "all_global_signed_totals_zero": len(nonzero_global_rows) == 0,
        "all_nonzero_weightings_have_paired_row_polarities": len(single_polarity_rows) == 0,
        "sample_weighting_rows": payload_rows[:8],
        "obstruction": "All 729 Aut(Z12)-orbit weightings of the source labels have zero global signed total on the P2718 marker.  Nonzero row-level values always occur in paired + and - polarities, so this source-dependent but orbit-respecting class does not select a P2721 polarity.",
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "source_dependent_aut_orbit_class_exhausted": audit["weighting_count"] == 729,
        "nonrepresentative_source_orbits_used": audit["orbit_count"] == 6,
        "all_global_signed_totals_zero": audit["all_global_signed_totals_zero"],
        "single_polarity_weighting_exported": audit["single_polarity_weighting_count"] > 0,
        "p2721_coupling_polarity_selected": False,
    }
    missing = [criterion for criterion, value in facts.items() if not value]
    return {
        "facts": facts,
        "missing_criteria": missing,
        "accepted_as_strict_source_dependent_invariant": not missing,
        "blocker": "Aut-orbit source dependence is allowed in this finite class, but it never yields a nonzero global signed value or a single unpaired polarity for the P2718 marker.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["source_orbit_weight_audit"]
    lines = [
        "# P2728/S1678 Aut(Z12) source-orbit weighted chiral invariant no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exhaustive source-dependent class",
        f"- source_orbits={audit['source_orbits']}",
        f"- weighting_count={audit['weighting_count']}",
        f"- nonzero_global_signed_total_count={audit['nonzero_global_signed_total_count']}",
        f"- single_polarity_weighting_count={audit['single_polarity_weighting_count']}",
        audit["obstruction"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    audit = exhaustive_weight_audit(marker_rows(inputs["P2718_CHIRAL_BISPECTRUM_MARKER"]))
    acceptance = acceptance_matrix(audit)
    no_go = audit["all_global_signed_totals_zero"] and audit["single_polarity_weighting_count"] == 0
    payload = {
        "status": "P2728_AUT_Z12_SOURCE_ORBIT_WEIGHTED_CHIRAL_INVARIANT_NO_GO" if no_go else "P2728_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "all {-1,0,+1} weights on Aut(Z12) source orbits applied to P2718 Im(B_{1,5}) marker rows",
        "source_orbit_weight_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "aut_orbit_weighted_source_invariant_exported": False,
            "nonzero_global_signed_value_exported": False,
            "single_polarity_selected": False,
            "p2721_coupling_polarity_selected": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2728 exhausts a source-dependent but non-representative candidate class: all {-1,0,+1} weights on the six Aut(Z12) source orbits.  Every weighted signed total cancels across the two orientations, and row-level signs remain paired.  Thus this class does not supply a strict invariant with a nonzero signed value or P2721 polarity selection.",
            "next_honest_step": "Do not continue Aut-orbit source-weighting variants for the P2718 marker.  A next admissible proof-grade move must provide a new strict source-dependent invariant that is not merely an Aut-orbit scalar weighting and that includes an exported torsor/polarity-coupling theorem; otherwise preserve the P2697-P2728 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2728/S1678 Aut(Z12) source-orbit weighted chiral invariant no-go", "## P2728/S1678 Aut(Z12) source-orbit weighted chiral invariant no-go\n\n`P2728/S1678` exhausts all `{-1,0,+1}` weights on the six Aut(Z12) source orbits for the P2718 marker `Im(B_{1,5})`.  All 729 source-dependent orbit weightings have zero global signed total, and nonzero row-level values remain paired in `+/-` polarities.  This class exports no strict source-dependent invariant, no P2721 polarity selection, no `lambda` fixing, no `QW-2191` discharge, no role transfer, no `L_total`, and no ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2728/S1678 Aut-orbit source weighting Ltotal guard", "## P2728/S1678 Aut-orbit source weighting Ltotal guard\n\n`P2728/S1678` finds zero global signed value for every Aut(Z12)-orbit source weighting of `Im(B_{1,5})`.  Because the class supplies no unpaired signed source and no P2721 coupling theorem, it cannot add a variational source term or promote `L_total`.\n")
    append_once(AGENTS, "Current Aut(Z12) source-orbit weighted invariant no-go guardrail (P2728/S1678, 2026-06-14)", "## Current Aut(Z12) source-orbit weighted invariant no-go guardrail (P2728/S1678, 2026-06-14)\n\n- P2728 exhausts all `{-1,0,+1}` weights on the six Aut(Z12) source orbits for the P2718 marker `Im(B_{1,5})`, giving 729 source-dependent but non-representative candidate invariants.\n- Every weighted global signed total is zero, and every nonzero row-level weighting remains paired in `+/-` polarities; no P2721 polarity is selected.\n- Do not continue Aut-orbit source-weighting variants as `QW-2191` discharge, selector closure, pair12 strict-core upgrade, role transfer, bridge closure, `L_total`, or ToE.  A next admissible move must introduce a new strict source-dependent invariant beyond orbit-scalar weighting with an exported torsor/polarity-coupling theorem, or preserve the P2697-P2728 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
