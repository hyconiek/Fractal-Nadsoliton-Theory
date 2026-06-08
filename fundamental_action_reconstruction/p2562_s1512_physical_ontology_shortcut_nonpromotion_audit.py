#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2562_s1512_physical_ontology_shortcut_nonpromotion_audit.json"
MD = GEN / "p2562_s1512_physical_ontology_shortcut_nonpromotion_audit.md"

SOURCE_FILES = {
    "P2561_POST_DAMPING_RESIDUAL_BRIDGE": GEN / "p2561_s1511_strict_completion_post_damping_residual_bridge_two_key_certificate.json",
    "P2515_OPERATOR_ORDER_SIGNATURE": GEN / "p2515_s1465_strict_damping_rg_operator_order_signature_acceptance_audit.json",
    "N50_KERNEL_NONIDENTIFICATION": ROOT / "n50_current_legacy_ontological_kernel_to_strict_gate_kernel_nonidentification_theorem.py",
    "F326_PHASE_FREQUENCY_OBSTRUCTION": ROOT / "F326_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_PHASE_FREQUENCY_NONCONFORMAL_OBSTRUCTION_WITNESS_PACKET.md",
}

STRICT = {
    "omega": Fraction(743, 4000),
    "phi": Fraction(13, 80),
    "eta": Fraction(9, 5),
}
LEGACY = {
    "omega": math.pi / 4.0,
    "phi": math.pi / 6.0,
    "linear_damping_power": Fraction(1, 1),
}
NEGATIVE_EXPORT_FLAGS = [
    "phenomenological_shortcut_promoted_to_source",
    "eta_fractal_tunneling_source_exported",
    "omega_nonlinear_resonance_source_exported",
    "hydrodynamic_m2_source_exported",
    "strict_damping_beta_eta_source_exported",
    "strict_dynamical_source_for_A_P_D_exported",
    "strict_phase_frequency_source_exported",
    "bridge_theorem_exported",
    "legacy_to_strict_completion_bridge_exported",
    "role_transfer_theorem_exported",
    "role_bearing_ltotal_exported",
    "qw2191_discharged_by_this_certificate",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "DIAGRAMS_KERNEL_TRANSFORMATION.md", "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2562|S1512|physical ontology shortcut|phenomenological shortcut|shortcut nonpromotion|ontology shortcut nonpromotion",
        "intended_research_nonduplication": "hydrodynamic.*m=2|m=2.*hydrodynamic|fractal tunneling.*eta|eta.*fractal tunneling|nonlinear resonance.*omega|omega.*nonlinear resonance",
        "legacy_phenomenology_inputs": "DIAGRAMS_KERNEL_TRANSFORMATION|K_total|fractal topology|hydrodynamic|lepkość|turbulent|winding|pi/4|pi/6",
        "phase_frequency_precursors": "F326|N50|phase/frequency nonconformal obstruction|explicit_phase_frequency_bridge_present|strict_phase_frequency_source",
        "operator_order_precursors": "P2515|S1465|operator-order signature|m=2|Euler-Lagrange order|strict_damping_beta_eta_source",
        "bridge_guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def file_presence(path: Any) -> dict[str, Any]:
    if hasattr(path, "exists") and path.exists():
        text = path.read_text(encoding="utf-8")
        return {"path": rel(path), "present": True, "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(), "bytes": len(text.encode("utf-8"))}
    return {"path": rel(path), "present": False}


def numerical_shortcut_audit() -> dict[str, Any]:
    strict_omega = float(STRICT["omega"])
    strict_phi = float(STRICT["phi"])
    strict_eta = float(STRICT["eta"])
    legacy_omega = LEGACY["omega"]
    legacy_phi = LEGACY["phi"]
    eta_gap = strict_eta - float(LEGACY["linear_damping_power"])
    rows = []
    for d in range(1, 12):
        legacy_phase = legacy_omega * d + legacy_phi
        strict_phase = strict_omega * d + strict_phi
        rows.append({
            "d": d,
            "legacy_phase": legacy_phase,
            "strict_phase": strict_phase,
            "phase_gap": legacy_phase - strict_phase,
            "legacy_cos": math.cos(legacy_phase),
            "strict_cos": math.cos(strict_phase),
            "cos_gap_abs": abs(math.cos(legacy_phase) - math.cos(strict_phase)),
        })
    integer_rescalings = []
    for k in range(1, 13):
        integer_rescalings.append({
            "integer_rescaling_k": k,
            "abs_legacy_omega_minus_k_strict_omega": abs(legacy_omega - k * strict_omega),
            "abs_legacy_phi_minus_k_strict_phi": abs(legacy_phi - k * strict_phi),
        })
    best_omega_integer_rescaling = min(integer_rescalings, key=lambda row: row["abs_legacy_omega_minus_k_strict_omega"])
    return {
        "strict_parameters": {"omega": str(STRICT["omega"]), "phi": str(STRICT["phi"]), "eta": str(STRICT["eta"])},
        "legacy_parameters": {"omega": "pi/4", "phi": "pi/6", "linear_damping_power": str(LEGACY["linear_damping_power"])},
        "eta_gap_strict_minus_legacy_linear_power": eta_gap,
        "omega_gap_abs": abs(legacy_omega - strict_omega),
        "phi_gap_abs": abs(legacy_phi - strict_phi),
        "phase_rows_d1_to_d11": rows,
        "max_phase_gap_abs_d1_to_d11": max(abs(row["phase_gap"]) for row in rows),
        "max_cos_gap_abs_d1_to_d11": max(row["cos_gap_abs"] for row in rows),
        "integer_rescaling_search_k1_to_k12": integer_rescalings,
        "best_omega_integer_rescaling_k1_to_k12": best_omega_integer_rescaling,
        "no_exact_legacy_to_strict_numeric_identity": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2561_payload = load_json(SOURCE_FILES["P2561_POST_DAMPING_RESIDUAL_BRIDGE"])
    p2515_payload = load_json(SOURCE_FILES["P2515_OPERATOR_ORDER_SIGNATURE"])
    p2561 = theorem(p2561_payload, "strict_completion_post_damping_residual_bridge_two_key_certificate")
    p2515 = theorem(p2515_payload, "strict_damping_rg_operator_order_signature_acceptance_audit")
    numeric = numerical_shortcut_audit()
    heuristic_claims = [
        {
            "claim": "legacy fractal/topological tunneling heuristically explains strict eta=9/5",
            "computational_check": "strict eta differs from the legacy linear damping power by 4/5",
            "source_status": "heuristic_only_no_strict_eta_source_export",
            "promotable_to_bridge_atom": False,
        },
        {
            "claim": "legacy resonance intuition heuristically explains strict omega=743/4000 and phi=13/80",
            "computational_check": "finite d=1..11 phase/cos gaps remain nonzero; N50/F326 still mark no explicit phase/frequency bridge",
            "source_status": "heuristic_only_no_strict_phase_frequency_source_export",
            "promotable_to_bridge_atom": False,
        },
        {
            "claim": "hydrodynamic/Laplacian analogy explains operator order m=2",
            "computational_check": "P2515 identifies the m=2 signature acceptance boundary but does not source it",
            "source_status": "heuristic_only_no_hydrodynamic_m2_source_export",
            "promotable_to_bridge_atom": False,
        },
    ]
    residual_atoms = p2561.get("residual_atoms_after_hypothetical_damping_source", [])
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2562_T1_physical_ontology_shortcut_nonpromotion_audit",
        "audited_chain": ["DIAGRAMS_KERNEL_TRANSFORMATION.md", "P2515/S1465", "P2561/S1511", "N50", "F326"],
        "p2561_residual_two_key_frontier_inherited": p2561.get("residual_atoms_are_jointly_required_after_damping") is True,
        "p2515_m2_signature_identified_not_sourced_inherited": p2515.get("p2506_roughness_m2_signature_identified_not_sourced") is True,
        "heuristic_physical_ontology_claims_audited": heuristic_claims,
        "heuristic_claim_count": len(heuristic_claims),
        "promotable_heuristic_claim_count": sum(1 for claim in heuristic_claims if claim["promotable_to_bridge_atom"]),
        "numerical_shortcut_audit": numeric,
        "residual_atoms_before_shortcut": residual_atoms,
        "residual_atoms_after_heuristic_shortcut": residual_atoms,
        "heuristic_shortcut_changes_bridge_truth_table": False,
        "physical_ontology_document_can_be_written_as_interpretation": True,
        "physical_ontology_document_can_replace_source_theorem": False,
        "recommended_next_honest_step": (
            "Convert one heuristic into an explicit source theorem before promotion: either derive the phase/frequency/topological-bit passage as a strict source under QW-2191, or derive the hydrodynamic m=2/operator-order principle with boundary data. Do not treat DIAGRAMS phenomenology as a bridge theorem by analogy alone."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2561_frontier_inherited": theorem_export["p2561_residual_two_key_frontier_inherited"],
        "p2515_m2_status_inherited": theorem_export["p2515_m2_signature_identified_not_sourced_inherited"],
        "three_heuristic_claims_audited": theorem_export["heuristic_claim_count"] == 3,
        "no_heuristic_claim_promoted": theorem_export["promotable_heuristic_claim_count"] == 0,
        "finite_phase_rows_are_d1_to_d11": len(numeric["phase_rows_d1_to_d11"]) == 11,
        "eta_gap_is_four_fifths": abs(numeric["eta_gap_strict_minus_legacy_linear_power"] - 0.8) < 1e-12,
        "shortcut_does_not_change_residual_atoms": theorem_export["residual_atoms_after_heuristic_shortcut"] == theorem_export["residual_atoms_before_shortcut"],
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    source_fingerprints = {
        "P2561_POST_DAMPING_RESIDUAL_BRIDGE": sha256_json(p2561_payload),
        "P2515_OPERATOR_ORDER_SIGNATURE": sha256_json(p2515_payload),
        "N50_KERNEL_NONIDENTIFICATION": file_presence(SOURCE_FILES["N50_KERNEL_NONIDENTIFICATION"]),
        "F326_PHASE_FREQUENCY_OBSTRUCTION": file_presence(SOURCE_FILES["F326_PHASE_FREQUENCY_OBSTRUCTION"]),
    }
    return {
        "packet_id": "P2562",
        "stage_id": "S1512",
        "status": "P2562_PHYSICAL_ONTOLOGY_SHORTCUT_NONPROMOTION_AUDIT_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "physical_ontology_shortcut_nonpromotion_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": source_fingerprints,
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["physical_ontology_shortcut_nonpromotion_audit"]["theorem_export"]
    n = t["numerical_shortcut_audit"]
    lines = [
        "# P2562/S1512 physical ontology shortcut nonpromotion audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- P2561 residual two-key frontier inherited: `{t['p2561_residual_two_key_frontier_inherited']}`.",
        f"- P2515 m=2 signature identified but not sourced inherited: `{t['p2515_m2_signature_identified_not_sourced_inherited']}`.",
        f"- Heuristic physical-ontology claims audited: `{t['heuristic_claim_count']}`.",
        f"- Heuristic claims promoted to bridge/source atoms: `{t['promotable_heuristic_claim_count']}`.",
        f"- Strict-minus-legacy damping-power gap: `{n['eta_gap_strict_minus_legacy_linear_power']}`.",
        f"- Absolute omega gap pi/4 vs 743/4000: `{n['omega_gap_abs']}`.",
        f"- Absolute phi gap pi/6 vs 13/80: `{n['phi_gap_abs']}`.",
        f"- Max finite phase gap on d=1..11: `{n['max_phase_gap_abs_d1_to_d11']}`.",
        f"- Max finite cosine gap on d=1..11: `{n['max_cos_gap_abs_d1_to_d11']}`.", "",
        "## Interpretation", "",
        "The DIAGRAMS phenomenology may be kept as interpretation, but this audit does not allow it to replace an explicit source theorem.  The residual bridge atoms remain unchanged after granting the physical-ontology shortcut as a heuristic.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No phenomenological shortcut, eta source, phase/frequency source, hydrodynamic m=2 source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "",
        "## Fingerprint", "", f"`{payload['physical_ontology_shortcut_nonpromotion_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2562/S1512` audits the tempting physical-ontology shortcut suggested by the legacy `DIAGRAMS_KERNEL_TRANSFORMATION.md` phenomenology: fractal tunneling for `eta=9/5`, nonlinear resonance for `omega=743/4000, phi=13/80`, and hydrodynamic/Laplacian intuition for `m=2`.  The finite computation records a strict-minus-legacy damping-power gap of `4/5` and nonzero phase/cosine gaps on `d=1..11`; it also inherits P2515's result that the `m=2` operator signature is identified but not sourced.  Therefore the shortcut may be documented as interpretation, but it does not export `strict_damping_beta_eta_source`, `strict_phase_frequency_source`, or a bridge theorem.
""".strip()
    lag_section = """
`P2562/S1512` blocks promotion of DIAGRAMS-style physical ontology directly into a role-bearing `L_total`.  Hydrodynamic and resonance analogies remain source targets, not source theorems; role-transfer, QW-2191, and ToE closure remain closed gates.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2562/S1512 physical ontology shortcut nonpromotion guard", "## P2562/S1512 physical ontology shortcut nonpromotion guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2562/S1512 physical ontology shortcut nonpromotion Ltotal guard", "## P2562/S1512 physical ontology shortcut nonpromotion Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
