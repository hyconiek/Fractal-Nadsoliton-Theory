#!/usr/bin/env python3
"""P2721/S1671: chiral-bispectrum sign to orientation-torsor coupling audit.

P2720 blocked source/phase-origin localization for the P2718 marker under the
Z12 translation orbit.  The other remaining P2719 obligation is a torsor
coupling theorem.  P2721 therefore asks a narrower finite question: does the
sign torsor carried by Im(B_{1,5}) admit an Aut(Z12)-equivariant coupling to
the P2708/P2714 orientation torsor, and does that coupling fix the missing
lambda sign without an additional polarity/source premise?
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.json"
MD = GEN / "p2721_s1671_chiral_bispectrum_sign_to_orientation_torsor_coupling_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

AUT_UNITS = (1, 5, 7, 11)
SIGN_TORSOR = (-2, 2)
ORIENTATION_TORSOR = ("-omega", "+omega")
INPUTS = {
    "P2720_TRANSLATION_ORBIT_LOCALIZER": GEN / "p2720_s1670_chiral_bispectrum_translation_orbit_phase_origin_localizer_no_go.json",
    "P2719_PHASE_ORIGIN_AUDIT": GEN / "p2719_s1669_chiral_bispectrum_phase_origin_source_localizer_theorem_audit.json",
    "P2718_CHIRAL_BISPECTRUM_FORMULA": GEN / "p2718_s1668_chiral_bispectrum_explicit_signed_formula_torsor_coupling_audit.json",
    "P2714_ORIENTATION_TORSOR": GEN / "p2714_s1664_z12_orientation_torsor_global_section_obstruction.json",
    "P2708_BOUNDARY_COCYCLE": GEN / "p2708_s1658_z12_boundary_1_cocycle_selector_source_obstruction.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "canonical_coupling_polarity_selected",
    "unconditional_torsor_coupling_theorem_exported",
    "phase_origin_source_localizer_exported",
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


def is_inversion_unit(unit: int) -> bool:
    return unit in (7, 11)


def act_on_sign(unit: int, sign: int) -> int:
    return -sign if is_inversion_unit(unit) else sign


def act_on_orientation(unit: int, orientation: str) -> str:
    if not is_inversion_unit(unit):
        return orientation
    return "+omega" if orientation == "-omega" else "-omega"


def enumerate_couplings() -> list[dict[str, Any]]:
    candidates = []
    for image_neg, image_pos in itertools.product(ORIENTATION_TORSOR, repeat=2):
        mapping = {-2: image_neg, 2: image_pos}
        failures = []
        for unit in AUT_UNITS:
            for sign in SIGN_TORSOR:
                lhs = mapping[act_on_sign(unit, sign)]
                rhs = act_on_orientation(unit, mapping[sign])
                if lhs != rhs:
                    failures.append({"unit": unit, "sign": sign, "lhs": lhs, "rhs": rhs})
        candidates.append({
            "mapping": {str(key): value for key, value in mapping.items()},
            "aut_equivariant": not failures,
            "failure_count": len(failures),
            "failures": failures,
            "polarity": "positive_marker_to_plus_omega" if mapping[2] == "+omega" else "positive_marker_to_minus_omega",
        })
    return candidates


def marker_sign_summary(p2718: dict[str, Any]) -> dict[str, Any]:
    rows = p2718.get("marker_rows", [])
    values_by_orientation: dict[str, set[float]] = {"-1": set(), "1": set()}
    for row in rows:
        values_by_orientation[str(row["orientation"])].add(float(row["marker_imag"]))
    return {
        "row_count": len(rows),
        "values_by_orientation": {key: sorted(value) for key, value in values_by_orientation.items()},
        "p2718_orientation_sign_pattern_recovered": values_by_orientation == {"-1": {2.0}, "1": {-2.0}},
    }


def coupling_audit(p2718: dict[str, Any], p2720: dict[str, Any]) -> dict[str, Any]:
    candidates = enumerate_couplings()
    equivariant = [candidate for candidate in candidates if candidate["aut_equivariant"]]
    localizer_exported = bool(p2720.get("translation_orbit_audit", {}).get("facts", {}).get("source_origin_localizer_exported"))
    polarity_selected = False
    facts = {
        "marker_sign_torsor_defined": marker_sign_summary(p2718)["p2718_orientation_sign_pattern_recovered"],
        "aut_equivariant_couplings_exist": len(equivariant) == 2,
        "coupling_polarity_pair_remains": len(equivariant) == 2,
        "phase_origin_source_localizer_exported": localizer_exported,
        "canonical_coupling_polarity_selected": polarity_selected,
        "unconditional_torsor_coupling_theorem_exported": localizer_exported and polarity_selected,
    }
    return {
        "marker_sign_summary": marker_sign_summary(p2718),
        "candidate_count": len(candidates),
        "aut_equivariant_coupling_count": len(equivariant),
        "candidate_couplings": candidates,
        "facts": facts,
        "accepted_unconditional_coupling": facts["unconditional_torsor_coupling_theorem_exported"],
        "blocker": "Two Aut-equivariant couplings from the chiral-bispectrum sign torsor to the orientation torsor exist, but they are exchanged by the unresolved coupling polarity.  Since P2720 exports no phase-origin/source localizer and no strict law selects the polarity, the coupling remains conditional and does not fix lambda.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2721/S1671 chiral-bispectrum sign to orientation-torsor coupling audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Coupling result",
        f"- aut_equivariant_coupling_count={payload['coupling_audit']['aut_equivariant_coupling_count']}",
        f"- accepted_unconditional_coupling={payload['coupling_audit']['accepted_unconditional_coupling']}",
        payload["coupling_audit"]["blocker"],
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
    audit = coupling_audit(inputs["P2718_CHIRAL_BISPECTRUM_FORMULA"], inputs["P2720_TRANSLATION_ORBIT_LOCALIZER"])
    no_unlock = not audit["accepted_unconditional_coupling"]
    payload = {
        "status": "P2721_CONDITIONAL_SIGN_TORSOR_COUPLINGS_EXIST_BUT_NO_CANONICAL_POLARITY" if no_unlock else "P2721_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_object": "Aut(Z12)-equivariant maps from the Im(B_{1,5}) sign torsor {-2,+2} to the P2708/P2714 orientation torsor {-omega,+omega}",
        "coupling_audit": audit,
        "decision": {
            "conditional_equivariant_coupling_exported": True,
            "canonical_coupling_polarity_selected": False,
            "unconditional_torsor_coupling_theorem_exported": False,
            "strict_chiral_bispectrum_source_exported": False,
            "strict_mechanism_fixing_lambda_exported": False,
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": "P2721 finds real representation-level progress: exactly two Aut(Z12)-equivariant couplings connect the P2718 chiral-bispectrum sign torsor to the P2708/P2714 orientation torsor.  However, they are the two opposite polarity choices.  Current artifacts still do not export a phase-origin/source localizer or a strict law selecting one polarity, so the result is conditional and does not fix lambda or discharge QW-2191.",
            "next_honest_step": "Do not promote the conditional coupling pair to closure.  A next admissible move must supply one strict polarity-selection/source law for this coupling, or a genuinely new origin-sensitive strict invariant; otherwise preserve the P2697-P2721 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2721/S1671 chiral-bispectrum sign to orientation-torsor coupling audit", "## P2721/S1671 chiral-bispectrum sign to orientation-torsor coupling audit\n\n`P2721/S1671` audits Aut(Z12)-equivariant couplings from the P2718 marker sign torsor `{-2,+2}` to the P2708/P2714 orientation torsor `{-omega,+omega}`.  Exactly two equivariant couplings exist, but they are opposite polarity choices.  Since P2720 exports no phase-origin/source localizer and current artifacts export no strict law selecting one polarity, the coupling remains conditional and does not fix `lambda`, discharge `QW-2191`, upgrade pair12 strict-core, close the bridge, start role transfer, promote `L_total`, or close ToE.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2721/S1671 sign-torsor coupling Ltotal guard", "## P2721/S1671 sign-torsor coupling Ltotal guard\n\n`P2721/S1671` identifies a conditional representation-level coupling pair, not a variational source construction.  Without a phase-origin/source localizer or a strict polarity-selection law, the chiral-bispectrum sign cannot promote `L_total`, selector closure, pair12 strict-core, role transfer, bridge closure, or ToE.\n")
    append_once(AGENTS, "Current chiral-bispectrum sign-torsor coupling guardrail (P2721/S1671, 2026-06-14)", "## Current chiral-bispectrum sign-torsor coupling guardrail (P2721/S1671, 2026-06-14)\n\n- P2721 audits Aut(Z12)-equivariant maps from the P2718 chiral-bispectrum sign torsor `{-2,+2}` to the P2708/P2714 orientation torsor `{-omega,+omega}`.\n- Exactly two equivariant couplings exist, but they are opposite polarity choices; current artifacts do not export a phase-origin/source localizer or a strict law selecting one polarity, so `lambda` remains unfixed and `QW-2191` remains open.\n- Do not promote the conditional coupling pair to selector closure, pair12 strict-core, role transfer, bridge closure, `L_total`, or ToE.  A next admissible move must supply one strict polarity-selection/source law for this coupling or a genuinely new origin-sensitive strict invariant; otherwise preserve the P2697-P2721 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
