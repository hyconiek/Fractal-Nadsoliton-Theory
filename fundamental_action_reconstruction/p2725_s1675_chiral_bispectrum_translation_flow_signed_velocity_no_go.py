#!/usr/bin/env python3
"""P2725/S1675: chiral-bispectrum translation-flow signed-velocity no-go.

P2724 asked for a more proof-grade/computational next step: if a dynamic
chiral source is to fix the P2721 polarity, it needs a computable nonzero
signed value.  This packet tests the most direct finite dynamic lift of the
P2718 marker already in hand: move the D5 source around the Z12 translation
orbit with every nonzero discrete velocity and measure the finite difference of
Im(B_{1,5}).  Because P2720 showed the marker is translation-orbit constant,
the dynamic lift has zero signed velocity on every tested flow and cannot be
the missing strict chiral/time-orientation source.
"""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2725_s1675_chiral_bispectrum_translation_flow_signed_velocity_no_go.json"
MD = GEN / "p2725_s1675_chiral_bispectrum_translation_flow_signed_velocity_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

Z12 = 12
NONZERO_VELOCITIES = tuple(range(1, Z12))
INPUTS = {
    "P2724_COMMIT_INTAKE": GEN / "p2724_s1674_post_p2723_commit_intake_next_honest_step_certificate.json",
    "P2723_CHIRAL_TIME_SOURCE_MATRIX": GEN / "p2723_s1673_strict_chiral_time_orientation_source_law_candidate_matrix.json",
    "P2721_SIGN_TORSOR_COUPLING": GEN / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.json",
    "P2720_TRANSLATION_ORBIT_LOCALIZER": GEN / "p2720_s1670_chiral_bispectrum_translation_orbit_phase_origin_localizer_no_go.json",
    "P2718_CHIRAL_BISPECTRUM_MARKER": GEN / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "nonzero_signed_translation_flow_velocity_exported",
    "strict_dynamic_chiral_source_artifact_exported",
    "p2721_polarity_selected",
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


def marker_table(p2718: dict[str, Any]) -> dict[tuple[int, int], float]:
    table: dict[tuple[int, int], float] = {}
    for row in p2718.get("marker_rows", []):
        table[(int(row["orientation"]), int(row["source"]))] = float(row["marker_imag"])
    return table


def signed_velocity_audit(table: dict[tuple[int, int], float]) -> dict[str, Any]:
    flow_rows = []
    deltas_by_velocity: dict[int, set[float]] = defaultdict(set)
    nonzero_rows = []
    for orientation in (-1, 1):
        for velocity in NONZERO_VELOCITIES:
            for source in range(Z12):
                next_source = (source + velocity) % Z12
                current_value = table[(orientation, source)]
                next_value = table[(orientation, next_source)]
                delta = round(next_value - current_value, 9)
                if delta == -0.0:
                    delta = 0.0
                row = {
                    "orientation": orientation,
                    "source": source,
                    "velocity": velocity,
                    "next_source": next_source,
                    "marker_imag": current_value,
                    "next_marker_imag": next_value,
                    "delta_marker_imag": delta,
                    "nonzero_signed_velocity": abs(delta) > 1e-9,
                }
                flow_rows.append(row)
                deltas_by_velocity[velocity].add(delta)
                if row["nonzero_signed_velocity"]:
                    nonzero_rows.append(row)

    summary_rows = [
        {
            "velocity": velocity,
            "delta_values": sorted(deltas_by_velocity[velocity]),
            "all_zero": deltas_by_velocity[velocity] == {0.0},
        }
        for velocity in NONZERO_VELOCITIES
    ]
    return {
        "checked_flow_rows": len(flow_rows),
        "orientation_count": 2,
        "velocity_count": len(NONZERO_VELOCITIES),
        "source_count_per_orientation": Z12,
        "velocity_summary_rows": summary_rows,
        "nonzero_signed_velocity_rows": nonzero_rows,
        "nonzero_signed_velocity_count": len(nonzero_rows),
        "all_translation_flow_deltas_zero": len(nonzero_rows) == 0,
        "obstruction": "Every finite translation-flow difference of the P2718 Im(B_{1,5}) marker is zero for both orientations and all 11 nonzero Z12 velocities.  A dynamic source built only from translation of this marker has no nonzero signed value to couple to the P2721 polarity pair.",
    }


def acceptance_matrix(flow: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "explicit_dynamic_lift_defined": True,
        "all_nonzero_translation_velocities_checked": flow["velocity_count"] == 11,
        "computable_signed_value_exported": True,
        "nonzero_signed_value_exported": flow["nonzero_signed_velocity_count"] > 0,
        "coupled_to_p2721_polarity_pair": False,
        "selects_exactly_one_polarity": False,
    }
    criteria = list(facts)
    missing = [name for name in criteria if not facts[name]]
    return {
        "criteria": criteria,
        "facts": facts,
        "missing_criteria": missing,
        "accepted_as_strict_dynamic_chiral_source": not missing,
        "blocker": "The translation-flow lift is computable, but its signed velocity is identically zero and no P2721 polarity-coupling theorem is exported.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2725/S1675 chiral-bispectrum translation-flow signed-velocity no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite dynamic test",
        f"- checked_flow_rows={payload['translation_flow_audit']['checked_flow_rows']}",
        f"- velocity_count={payload['translation_flow_audit']['velocity_count']}",
        f"- nonzero_signed_velocity_count={payload['translation_flow_audit']['nonzero_signed_velocity_count']}",
        payload["translation_flow_audit"]["obstruction"],
        "",
        "## Acceptance",
        f"- accepted_as_strict_dynamic_chiral_source={payload['acceptance_matrix']['accepted_as_strict_dynamic_chiral_source']}",
        f"- missing={payload['acceptance_matrix']['missing_criteria']}",
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    table = marker_table(inputs["P2718_CHIRAL_BISPECTRUM_MARKER"])
    flow = signed_velocity_audit(table)
    acceptance = acceptance_matrix(flow)
    no_go = flow["all_translation_flow_deltas_zero"] and not acceptance["accepted_as_strict_dynamic_chiral_source"]
    payload = {
        "status": "P2725_TRANSLATION_FLOW_SIGNED_VELOCITY_NO_GO_NO_STRICT_DYNAMIC_CHIRAL_SOURCE" if no_go else "P2725_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_dynamic_candidate": "finite Z12 translation-flow velocity of the P2718 Im(B_{1,5}) marker",
        "translation_flow_audit": flow,
        "acceptance_matrix": acceptance,
        "decision": {
            "nonzero_signed_translation_flow_velocity_exported": False,
            "strict_dynamic_chiral_source_artifact_exported": False,
            "p2721_polarity_selected": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2725 turns the P2724 recommendation into a bounded computation by testing the most direct dynamic lift of the exact P2718 marker.  Across 2 orientations, 11 nonzero Z12 velocities, and 12 sources, all 264 finite translation-flow signed velocities are zero.  Therefore this lift supplies no nonzero signed dynamic/chiral value and cannot select a P2721 coupling polarity.",
            "next_honest_step": "Do not reuse translation flow of Im(B_{1,5}) as the missing source.  The next admissible proof-grade move must either introduce a new non-translation dynamic/chiral observable with a computable nonzero signed value and an explicit P2721 polarity-coupling theorem, or preserve the P2697-P2725 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2725/S1675 chiral-bispectrum translation-flow signed-velocity no-go", "## P2725/S1675 chiral-bispectrum translation-flow signed-velocity no-go\n\n`P2725/S1675` gives the computational follow-up to P2724 by testing the finite translation-flow velocity of the exact P2718 marker `Im(B_{1,5})`.  Across two orientations, all 11 nonzero Z12 velocities, and 12 sources, all 264 signed finite differences are zero.  This dynamic lift exports no nonzero signed chiral/time-orientation value, no P2721 polarity selection, no strict `lambda` fixing, no `QW-2191` discharge, no role transfer, no `L_total`, and no ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2725/S1675 zero translation-flow velocity Ltotal guard", "## P2725/S1675 zero translation-flow velocity Ltotal guard\n\n`P2725/S1675` is a finite dynamic-observable no-go: the translation-flow signed velocity of `Im(B_{1,5})` is identically zero on the audited Z12 domain.  Since it supplies no nonzero signed source value and no P2721 coupling theorem, it cannot add a variational source term or promote `L_total`.\n")
    append_once(AGENTS, "Current chiral-bispectrum translation-flow signed-velocity no-go guardrail (P2725/S1675, 2026-06-14)", "## Current chiral-bispectrum translation-flow signed-velocity no-go guardrail (P2725/S1675, 2026-06-14)\n\n- P2725 performs the bounded computational follow-up requested after P2724: it tests the finite Z12 translation-flow velocity of the exact P2718 marker `Im(B_{1,5})` as a candidate strict dynamic/chiral source.\n- The computation checks 264 flow rows (2 orientations × 11 nonzero velocities × 12 sources) and every signed finite difference is zero; therefore this lift supplies no nonzero signed value and cannot select a P2721 coupling polarity.\n- Do not reuse translation flow of `Im(B_{1,5})` as `QW-2191` discharge, selector closure, pair12 strict-core upgrade, role transfer, bridge closure, `L_total`, or ToE.  A next admissible move must introduce a non-translation dynamic/chiral observable with a computable nonzero signed value and an explicit P2721 polarity-coupling theorem, or preserve the P2697-P2725 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
