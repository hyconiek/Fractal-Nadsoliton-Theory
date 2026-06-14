#!/usr/bin/env python3
"""P2719/S1669: phase-origin/source-localizer theorem audit for Im(B_{1,5}).

P2718 found a real signed chiral-bispectrum marker, but left exactly two
missing obligations: a non-premise phase-origin/source localizer and an exported
coupling theorem to the P2708/P2714 orientation torsor.  P2719 audits those
obligations for the exact formula rather than introducing another pseudoscalar
name or promoting the marker by wishful closure.
"""
from __future__ import annotations

import json
import hashlib
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2719_s1669_chiral_bispectrum_phase_origin_source_localizer_theorem_audit.json"
MD = GEN / "p2719_s1669_chiral_bispectrum_phase_origin_source_localizer_theorem_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2718_CHIRAL_BISPECTRUM_FORMULA": GEN / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.json",
    "P2717_SOURCE_MATRIX": GEN / "p2717_s1667_concrete_pseudoscalar_chiral_source_candidate_matrix.json",
    "P2716_PSEUDOSCALAR_ACCEPTANCE": GEN / "p2716_s1666_inversion_odd_pseudoscalar_source_acceptance_audit.json",
    "P2708_BOUNDARY_COCYCLE": GEN / "p2708_s1658_z12_boundary_1_cocycle_selector_source_obstruction.json",
    "P2714_ORIENTATION_TORSOR": GEN / "p2714_s1664_z12_orientation_torsor_global_section_obstruction.json",
    "P2367_PHASE_ORIGIN_NO_GO": GEN / "p2367_s1317_selector_phase_origin_admissibility_no_go_probe.json",
}

THEOREM_OBLIGATIONS = [
    "exact_formula_fixed",
    "signed_marker_nonzero_and_orientation_separating",
    "phase_origin_reference_internal_to_strict_artifacts",
    "phase_origin_not_translation_gauge_or_source_label",
    "source_localizer_selects_one_origin_nonpremise",
    "torsor_coupling_theorem_exported",
    "coupling_is_qw2191_safe_and_nonselector_replay",
]

NEGATIVE_EXPORT_FLAGS = [
    "phase_origin_source_localizer_theorem_exported",
    "strict_chiral_bispectrum_source_exported",
    "torsor_coupling_theorem_exported",
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


def audit_phase_origin(p2718: dict[str, Any]) -> dict[str, Any]:
    summary = p2718.get("finite_summary", {})
    acceptance = p2718.get("acceptance_matrix", {}).get("facts", {})
    exact_formula_fixed = "Im(B_{1,5})" in summary.get("formula", "")
    signed_ok = bool(summary.get("orientation_separating")) and bool(acceptance.get("nonzero_signed_value_on_all_rows"))
    source_blind = bool(summary.get("source_blind_under_translation"))

    facts = {
        "exact_formula_fixed": exact_formula_fixed,
        "signed_marker_nonzero_and_orientation_separating": signed_ok,
        "phase_origin_reference_internal_to_strict_artifacts": False,
        "phase_origin_not_translation_gauge_or_source_label": False,
        "source_localizer_selects_one_origin_nonpremise": False,
        "torsor_coupling_theorem_exported": False,
        "coupling_is_qw2191_safe_and_nonselector_replay": False,
    }
    missing = [name for name in THEOREM_OBLIGATIONS if not facts[name]]
    rows = [
        {
            "obligation": "exact chiral-bispectrum formula",
            "result": exact_formula_fixed,
            "comment": "P2718 fixes the audited formula as Im(B_{1,5}); P2719 does not introduce a new pseudoscalar candidate.",
        },
        {
            "obligation": "nonzero signed marker",
            "result": signed_ok,
            "comment": "P2718 supplies the finite +2/-2 orientation-separating witness on all 24 rows.",
        },
        {
            "obligation": "non-premise phase-origin/source localizer",
            "result": False,
            "comment": "Current artifacts still treat the origin/source index as translation/source-label data; no internal strict theorem selects one phase origin non-premise.",
        },
        {
            "obligation": "torsor coupling theorem",
            "result": False,
            "comment": "No exported theorem maps the chiral-bispectrum sign to the P2708/P2714 orientation torsor as a QW-2191-safe strict source.",
        },
    ]
    return {
        "facts": facts,
        "missing_obligations": missing,
        "source_blind_under_translation_used": source_blind,
        "phase_origin_source_localizer_theorem_exported": not missing,
        "audit_rows": rows,
        "blocker": "The formula-level sign survives, but the origin of the phase remains source-label/translation-gauge rather than a non-premise strict localizer, and the torsor-coupling theorem is still absent.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2719/S1669 chiral-bispectrum phase-origin/source-localizer theorem audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        f"- theorem_exported={payload['theorem_audit']['phase_origin_source_localizer_theorem_exported']}",
        f"- missing={payload['theorem_audit']['missing_obligations']}",
        payload["theorem_audit"]["blocker"],
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
    theorem_audit = audit_phase_origin(inputs["P2718_CHIRAL_BISPECTRUM_FORMULA"])
    no_unlock = not theorem_audit["phase_origin_source_localizer_theorem_exported"]
    payload = {
        "status": "P2719_CHIRAL_BISPECTRUM_PHASE_ORIGIN_LOCALIZER_THEOREM_NO_UNLOCK" if no_unlock else "P2719_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_formula": "Im(B_{1,5}) chiral-bispectrum marker from P2718",
        "theorem_obligations": THEOREM_OBLIGATIONS,
        "theorem_audit": theorem_audit,
        "decision": {
            "phase_origin_source_localizer_theorem_exported": False,
            "strict_chiral_bispectrum_source_exported": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2719 audits the exact P2718 chiral-bispectrum marker.  The finite signed marker remains real, but the repository still lacks a non-premise phase-origin/source localizer and lacks an exported theorem coupling that sign to the P2708/P2714 orientation torsor.  Therefore the marker is not promoted to a strict source fixing lambda or discharging QW-2191.",
            "next_honest_step": "Do not replay generic pseudoscalar enumeration or promote Im(B_{1,5}) by itself.  A further move is admissible only if it supplies one explicit non-premise phase-origin/source-localizer theorem or one exported torsor-coupling theorem for this exact formula; otherwise preserve the P2697-P2719 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2719/S1669 chiral-bispectrum phase-origin theorem audit", "## P2719/S1669 chiral-bispectrum phase-origin theorem audit\n\n`P2719/S1669` audits the exact P2718 marker `Im(B_{1,5})` for the missing phase-origin/source-localizer and torsor-coupling theorem.  The finite signed formula evidence remains intact, but current artifacts do not export a non-premise phase-origin/source localizer and do not export a theorem coupling the sign to the P2708/P2714 orientation torsor.  No strict `lambda` fixing, `QW-2191` discharge, pair12 strict-core upgrade, bridge closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2719/S1669 chiral-bispectrum phase-origin Ltotal guard", "## P2719/S1669 chiral-bispectrum phase-origin Ltotal guard\n\n`P2719/S1669` is a theorem-obligation audit, not a variational source construction.  Since the non-premise phase-origin/source localizer and P2708/P2714 torsor-coupling theorem remain absent, the chiral-bispectrum sign cannot promote `L_total`, selector closure, pair12 strict-core, role transfer, bridge closure, or ToE.\n")
    append_once(AGENTS, "Current chiral-bispectrum phase-origin/source-localizer theorem guardrail (P2719/S1669, 2026-06-14)", "## Current chiral-bispectrum phase-origin/source-localizer theorem guardrail (P2719/S1669, 2026-06-14)\n\n- P2719 audits the exact P2718 marker `Im(B_{1,5})` for the missing non-premise phase-origin/source localizer and the missing coupling theorem to the P2708/P2714 orientation torsor.\n- The finite signed marker remains real, but current artifacts do not export the required localizer or torsor-coupling theorem; therefore the marker still does not fix `lambda` or discharge `QW-2191`.\n- Do not replay generic pseudoscalar enumeration or promote the chiral-bispectrum marker alone to selector closure, pair12 strict-core, role transfer, bridge closure, `L_total`, or ToE.  A next admissible move must supply one explicit non-premise phase-origin/source-localizer theorem or one exported torsor-coupling theorem for this exact formula; otherwise preserve the P2697-P2719 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
