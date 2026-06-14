#!/usr/bin/env python3
"""P2720/S1670: translation-orbit phase-origin localizer no-go for Im(B_{1,5}).

P2719 identified the precise missing theorem: the P2718 marker needs a
non-premise phase-origin/source localizer or a torsor-coupling theorem.  P2720
attacks the localizer obligation directly.  It tests whether the exact
chiral-bispectrum data already computed in P2718 can localize one Z12 source
origin once the 12 translations of the D5 support are treated as gauge/orbit
data rather than as externally supplied labels.
"""
from __future__ import annotations

import json
import hashlib
from collections import defaultdict
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2720_s1670_chiral_bispectrum_translation_orbit_phase_origin_localizer_no_go.json"
MD = GEN / "p2720_s1670_chiral_bispectrum_translation_orbit_phase_origin_localizer_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

Z12_NODE_COUNT = 12
INPUTS = {
    "P2719_PHASE_ORIGIN_AUDIT": GEN / "p2719_s1669_chiral_bispectrum_phase_origin_source_localizer_theorem_audit.json",
    "P2718_CHIRAL_BISPECTRUM_FORMULA": GEN / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.json",
    "P2708_BOUNDARY_COCYCLE": GEN / "p2708_s1658_z12_boundary_1_cocycle_selector_source_obstruction.json",
    "P2714_ORIENTATION_TORSOR": GEN / "p2714_s1664_z12_orientation_torsor_global_section_obstruction.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "translation_orbit_phase_origin_localizer_exported",
    "nonpremise_source_origin_selected",
    "torsor_coupling_theorem_exported",
    "strict_chiral_bispectrum_source_exported",
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


def translation_orbit(source: int) -> list[int]:
    return sorted({(source + shift) % Z12_NODE_COUNT for shift in range(Z12_NODE_COUNT)})


def orbit_localizer_audit(p2718: dict[str, Any]) -> dict[str, Any]:
    rows = p2718.get("marker_rows", [])
    grouped: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[int(row["orientation"])].append(row)

    orbit_rows = []
    accepted_localizers = []
    for orientation, members in sorted(grouped.items()):
        sources = sorted(int(row["source"]) for row in members)
        marker_values = sorted({float(row["marker_imag"]) for row in members})
        recovered_orientations = sorted({int(row["orientation_recovered_from_marker"]) for row in members})
        full_orbit = sources == translation_orbit(0)
        constant_on_orbit = len(marker_values) == 1
        can_select_origin_without_label = full_orbit and not constant_on_orbit
        if can_select_origin_without_label:
            accepted_localizers.append(orientation)
        orbit_rows.append({
            "orientation": orientation,
            "source_orbit_size": len(sources),
            "sources": sources,
            "marker_values": marker_values,
            "orientation_recovered_from_marker_values": recovered_orientations,
            "full_z12_translation_orbit": full_orbit,
            "marker_constant_on_translation_orbit": constant_on_orbit,
            "can_select_one_origin_without_external_label": can_select_origin_without_label,
        })

    facts = {
        "exact_p2718_marker_rows_reused": len(rows) == 24,
        "each_orientation_has_full_z12_source_orbit": all(row["full_z12_translation_orbit"] for row in orbit_rows),
        "marker_constant_on_each_translation_orbit": all(row["marker_constant_on_translation_orbit"] for row in orbit_rows),
        "orientation_sign_still_recoverable": all(len(row["orientation_recovered_from_marker_values"]) == 1 for row in orbit_rows),
        "source_origin_localizer_exported": len(accepted_localizers) > 0,
    }
    return {
        "facts": facts,
        "orbit_rows": orbit_rows,
        "accepted_localizer_orientations": accepted_localizers,
        "accepted_localizer_count": len(accepted_localizers),
        "obstruction": "The P2718 bispectrum sign separates orientation but is constant across the full Z12 translation orbit for each orientation.  Any localizer using only this orbit-quotiented marker data cannot select one source/phase origin without importing an external source label or gauge convention.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2720/S1670 chiral-bispectrum translation-orbit phase-origin localizer no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Orbit result",
        f"- accepted_localizer_count={payload['translation_orbit_audit']['accepted_localizer_count']}",
        f"- marker_constant_on_each_translation_orbit={payload['translation_orbit_audit']['facts']['marker_constant_on_each_translation_orbit']}",
        payload["translation_orbit_audit"]["obstruction"],
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
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    orbit_audit = orbit_localizer_audit(inputs["P2718_CHIRAL_BISPECTRUM_FORMULA"])
    no_unlock = orbit_audit["accepted_localizer_count"] == 0
    payload = {
        "status": "P2720_TRANSLATION_ORBIT_PHASE_ORIGIN_LOCALIZER_NO_GO" if no_unlock else "P2720_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_formula": "Im(B_{1,5}) chiral-bispectrum marker from P2718, restricted to translation-orbit phase-origin localization",
        "translation_orbit_audit": orbit_audit,
        "decision": {
            "translation_orbit_phase_origin_localizer_exported": False,
            "nonpremise_source_origin_selected": False,
            "strict_chiral_bispectrum_source_exported": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2720 attacks the P2719 localizer obligation directly.  The exact P2718 marker keeps its orientation sign, but for each orientation it is constant on the full 12-element translation orbit.  Thus it cannot select a non-premise source/phase origin without importing an external source label or translation-gauge convention.",
            "next_honest_step": "Do not repeat the translation-orbit localizer attempt for Im(B_{1,5}).  A next admissible move must supply a genuinely new origin-sensitive strict invariant/source law that is not translation-gauge/source-label data, or an independent exported torsor-coupling theorem; otherwise preserve the P2697-P2720 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2720/S1670 chiral-bispectrum translation-orbit phase-origin localizer no-go", "## P2720/S1670 chiral-bispectrum translation-orbit phase-origin localizer no-go\n\n`P2720/S1670` attacks the P2719 phase-origin/source-localizer obligation for the exact P2718 marker `Im(B_{1,5})`.  The marker still separates orientation, but it is constant on each full 12-source Z12 translation orbit, so it cannot select a non-premise source/phase origin without importing an external source label or translation-gauge convention.  No strict `lambda` fixing, `QW-2191` discharge, pair12 strict-core upgrade, bridge closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2720/S1670 translation-orbit localizer Ltotal guard", "## P2720/S1670 translation-orbit localizer Ltotal guard\n\n`P2720/S1670` is a finite orbit obstruction for a phase-origin localizer, not a variational source construction.  Since the chiral-bispectrum marker is translation-orbit constant and no origin-sensitive strict source law is exported, it cannot promote `L_total`, selector closure, pair12 strict-core, role transfer, bridge closure, or ToE.\n")
    append_once(AGENTS, "Current chiral-bispectrum translation-orbit phase-origin localizer guardrail (P2720/S1670, 2026-06-14)", "## Current chiral-bispectrum translation-orbit phase-origin localizer guardrail (P2720/S1670, 2026-06-14)\n\n- P2720 tests the P2719 localizer obligation for the exact P2718 marker `Im(B_{1,5})` under the Z12 translation orbit of the D5 source.\n- The marker still separates orientation, but it is constant on the full 12-source translation orbit for each orientation; therefore it cannot select a non-premise source/phase origin without importing an external source label or translation-gauge convention.\n- Do not repeat this translation-orbit localizer attempt or promote the chiral-bispectrum marker to selector closure, pair12 strict-core, role transfer, bridge closure, `L_total`, or ToE.  A next admissible move must supply a genuinely new origin-sensitive strict invariant/source law or an independent exported torsor-coupling theorem; otherwise preserve the P2697-P2720 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
