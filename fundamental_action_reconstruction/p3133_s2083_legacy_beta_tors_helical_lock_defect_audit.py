"""P3133/S2083: legacy beta_tors helical-lock defect audit.

P3132 sharpened the only admissible reopening condition for the interlocked
helix lane: exhibit one strict formula/artifact for a nontranslation,
support-local helical-lock defect D_HL with a nonzero inversion-odd value and
coupling to one support representative.

This audit tests the user's legacy-kernel suspicion directly: did the original
sinusoid plus beta_tors denominator already contain that missing D_HL, or was
part of the torsion/phase semantics lost in the strict operational refreeze?
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p3133_s2083_legacy_beta_tors_helical_lock_defect_audit.json"
MD = GEN / "p3133_s2083_legacy_beta_tors_helical_lock_defect_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
P3132 = GEN / "p3132_s2082_interlocked_helix_support_local_section_audit.json"
K1 = ROOT / "K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md"
K2 = ROOT / "K2_STRICT_GATE_KERNEL_DERIVATION_CHAIN_NOTE.md"
F2 = ROOT / "F2_STRICT_GATE_KERNEL_PROVENANCE_AND_FAR_INPUT_CLASSIFICATION_PACKET.md"
S2 = ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"
QW1729 = REPO / "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW1729_NADSOLITON_KERNEL_CHARACTERISTICS_MAP.md"
QW2041 = REPO / "material_dowodowy/korpus_qw_pozostaly/raporty_md/RAPORT_QW2041_CANONICAL_REFROZEN_REPARAMETERIZATION_AUDIT.md"


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text(encoding="utf-8") if path.exists() else ""
    if marker not in existing:
        path.write_text(existing.rstrip() + "\n\n" + text.strip() + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    rows = [
        {
            "legacy_atom": "alpha_geo amplitude",
            "legacy_location": "multiplicative numerator normalization",
            "strict_successor": "implicit/absorbed operational normalization in K_strict_gate",
            "D_HL_nontranslation_defect": False,
            "nonzero_inversion_odd_value": False,
            "support_representative_coupling": False,
            "fresh_gap": "amplitude can scale a kernel but cannot by itself choose a helical support origin or chirality polarity",
        },
        {
            "legacy_atom": "cos(omega*d+phi) sinusoid",
            "legacy_location": "oscillatory phase/resonance carrier",
            "strict_successor": "omega,phi refrozen by the QW2038-QW2050 spectral/intersection chain",
            "D_HL_nontranslation_defect": True,
            "nonzero_inversion_odd_value": False,
            "support_representative_coupling": False,
            "fresh_gap": "the sinusoid preserves a phase/winding trace, but cos is even and translation-orbit blind unless a separate phase-origin localizer is exported",
        },
        {
            "legacy_atom": "beta_tors linear denominator",
            "legacy_location": "1/(1+beta_tors*d) torsion/topological path-summation damping",
            "strict_successor": "1/(1+beta*d**eta) nonlinear damping/compression with beta=1, eta=1.8 in the working strict gate",
            "D_HL_nontranslation_defect": False,
            "nonzero_inversion_odd_value": False,
            "support_representative_coupling": False,
            "fresh_gap": "beta_tors is a scalar damping/torsion marker; current artifacts do not make it an inversion-odd, support-local helical-lock defect",
        },
        {
            "legacy_atom": "combined alpha_geo*cos/(1+beta_tors*d)",
            "legacy_location": "effective ontological bridge kernel",
            "strict_successor": "completed/enriched strict gate only under explicit completion-map evidence",
            "D_HL_nontranslation_defect": True,
            "nonzero_inversion_odd_value": False,
            "support_representative_coupling": False,
            "fresh_gap": "the combined formula can motivate D_HL, but it still lacks an exported odd sign value and an absolute support representative after the diagonal Z12 quotient",
        },
    ]
    gates = [
        "strict nadsoliton provenance for D_HL rather than legacy-role transfer",
        "nontranslation defect not erased by common phase translation",
        "nonzero inversion-odd signed value",
        "coupling to one support representative after quotienting",
        "torsor polarity law rather than a chosen convention",
        "no selector/apparatus/observed-light/Planck/Lagrangian-normalization/bridge-role-transfer import",
    ]
    gate_matrix = []
    for row in rows:
        gate_matrix.append({
            "legacy_atom": row["legacy_atom"],
            "passed_gates": sum([
                row["D_HL_nontranslation_defect"],
                row["nonzero_inversion_odd_value"],
                row["support_representative_coupling"],
            ]),
            "required_core_gates": 3,
            "accepted_D_HL_source": row["D_HL_nontranslation_defect"] and row["nonzero_inversion_odd_value"] and row["support_representative_coupling"],
        })
    return {
        "status": "P3133_LEGACY_BETA_TORS_HELICAL_LOCK_DEFECT_GAP_CERTIFICATE",
        "audit_question": "Does the legacy beta_tors plus sinusoid kernel already export the P3132-required strict helical-lock defect D_HL?",
        "legacy_kernel": "K_legacy_ont(d)=alpha_geo*cos(omega*d+phi)/(1+beta_tors*d)",
        "strict_kernel": "K_strict_gate(d)=cos(omega*d+phi)/(1+beta*d^eta)",
        "input_hashes": {"P3132": sha(P3132), "K1": sha(K1), "K2": sha(K2), "F2": sha(F2), "S2": sha(S2), "QW1729": sha(QW1729), "QW2041": sha(QW2041)},
        "source_availability_note": "The external DIAGRAMS_KERNEL_TRANSFORMATION.md path was not available in this container; the audit therefore uses repo-restored legacy-kernel packets and QW1729/QW2041 provenance records.",
        "comparison_rows": rows,
        "acceptance_gates": gates,
        "gate_matrix": gate_matrix,
        "finite_certificate": {
            "legacy_atoms_tested": len(rows),
            "accepted_D_HL_sources": sum(row["accepted_D_HL_source"] for row in gate_matrix),
            "sinusoid_retains_phase_trace": True,
            "beta_tors_retains_scalar_torsion_damping_trace": True,
            "legacy_role_transfer_started": False,
            "strict_bridge_completion_exported": False,
        },
        "decision": {
            "bounded_result": "The user's suspicion is partly right: something important became under-visible in the move from the legacy formula to the strict operational gate. The legacy numerator keeps an explicit sinusoidal phase/resonance carrier and beta_tors keeps an explicit torsion-damping label, whereas the strict refreeze records omega, phi, beta, eta as operational gate parameters. However, this is not yet the missing P3132 D_HL source. The legacy formula provides motivation and a gap diagnosis, but current artifacts do not export beta_tors or the cosine phase as a nonzero inversion-odd, support-local helical-lock defect coupled to an absolute support representative.",
            "positive_scoped_flags": {
                "legacy_phase_trace_identified": True,
                "legacy_beta_tors_torsion_trace_identified": True,
                "strict_operational_absorption_gap_identified": True,
                "D_HL_acceptance_obligations_sharpened": True,
            },
            "negative_export_flags": {
                "D_HL_source_exported": False,
                "Zeta_OS_exported": False,
                "Gamma_SO_exported": False,
                "QW_2191_discharged": False,
                "legacy_role_transfer_exported": False,
                "bridge_completion_exported": False,
                "L_total_exported": False,
                "ToE_closure_exported": False,
            },
            "next_honest_step": "Build exactly one explicit candidate formula for D_HL from the legacy phase/torsion split, for example an oriented phase-gradient/torsion residual that is odd under inversion and is evaluated on a support-local representative. Then audit only its provenance, diagonal-translation behavior, sign polarity, and support coupling before any Zeta_OS/Gamma_SO retest. Do not start physical role transfer or claim bridge completion from beta_tors alone.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["finite_certificate"]
    decision = payload["decision"]
    lines = [
        "# P3133/S2083 Legacy beta_tors helical-lock defect audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Question",
        payload["audit_question"],
        "",
        "## Fresh comparison",
        "",
    ]
    for row in payload["comparison_rows"]:
        lines.extend([
            f"### {row['legacy_atom']}",
            f"- legacy role: {row['legacy_location']}",
            f"- strict successor: {row['strict_successor']}",
            f"- D_HL-like nontranslation trace: `{row['D_HL_nontranslation_defect']}`",
            f"- nonzero inversion-odd value: `{row['nonzero_inversion_odd_value']}`",
            f"- support-representative coupling: `{row['support_representative_coupling']}`",
            f"- gap: {row['fresh_gap']}",
            "",
        ])
    lines.extend([
        "## Certificate",
        f"- legacy atoms tested: `{cert['legacy_atoms_tested']}`",
        f"- accepted D_HL sources: `{cert['accepted_D_HL_sources']}`",
        f"- sinusoid retains phase trace: `{cert['sinusoid_retains_phase_trace']}`",
        f"- beta_tors retains scalar torsion/damping trace: `{cert['beta_tors_retains_scalar_torsion_damping_trace']}`",
        f"- bridge completion exported: `{cert['strict_bridge_completion_exported']}`",
        "",
        "## Decision",
        decision["bounded_result"],
        "",
        "## Recommendation",
        decision["next_honest_step"],
        "",
        "## Source availability note",
        payload["source_availability_note"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3133/S2083 legacy beta_tors helical-lock defect audit", "## P3133/S2083 legacy beta_tors helical-lock defect audit\n\n`P3133/S2083` executes the P3132-requested `D_HL` intake from the legacy-kernel side. It compares `alpha_geo`, `cos(omega*d+phi)`, `beta_tors`, and the combined legacy kernel against the strict helical-lock defect gates. The audit finds a real under-visibility gap: the legacy sinusoid preserves an explicit phase/resonance trace and `beta_tors` preserves an explicit torsion-damping label, while the strict gate stores these as operational parameters. But `0/4` tested legacy atoms export an accepted `D_HL` source because none supplies both a nonzero inversion-odd value and a support-representative coupling after the diagonal `Z12` quotient. No `Zeta_OS`, `Gamma_SO`, selector closure, bridge completion, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3133/S2083 legacy torsion is not yet a variational D_HL source", "## P3133/S2083 legacy torsion is not yet a variational D_HL source\n\n`P3133/S2083` identifies a legacy phase/torsion under-visibility gap, but it does not export a Lagrangian density, Hamiltonian normalization, spacetime EOM, physical unit, `L_total`, bridge-completion theorem, or role-transfer theorem.\n")
    append_once(AGENTS, "Current legacy beta_tors helical-lock defect guardrail (P3133/S2083, 2026-07-12)", "## Current legacy beta_tors helical-lock defect guardrail (P3133/S2083, 2026-07-12)\n\n- P3133 tests whether the restored legacy kernel ingredients, especially `cos(omega*d+phi)` and `beta_tors`, already export the P3132-required strict helical-lock defect `D_HL`.\n- The audit confirms a real legacy-to-strict under-visibility gap: the legacy sinusoid carries explicit phase/resonance structure and `beta_tors` carries explicit torsion-damping semantics, while the strict operational gate absorbs them into `omega/phi/beta/eta`; nevertheless `0/4` tested legacy atoms export an accepted `D_HL` source because none provides both a nonzero inversion-odd value and a support-representative coupling after the diagonal `Z12` quotient.\n- Do not promote `beta_tors`, the cosine phase, or the combined legacy kernel to `Zeta_OS`, `Gamma_SO`, `QW-2191` discharge, bridge completion, role transfer, `L_total`, or ToE closure.\n- The next admissible move is exactly one explicit candidate formula for `D_HL` built from the legacy phase/torsion split, audited for provenance, diagonal-translation behavior, sign polarity, and support coupling before any retest.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
