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
OUT = GEN / "p2564_s1514_phase_frequency_finite_sign_cell_nonidentifiability_certificate.json"
MD = GEN / "p2564_s1514_phase_frequency_finite_sign_cell_nonidentifiability_certificate.md"

SOURCE_FILES = {
    "P2563_RATIONAL_WINDING_OBSTRUCTION": GEN / "p2563_s1513_phase_frequency_rational_winding_quotient_obstruction_certificate.json",
    "P2561_POST_DAMPING_RESIDUAL_BRIDGE": GEN / "p2561_s1511_strict_completion_post_damping_residual_bridge_two_key_certificate.json",
    "P2415_PHASE_AFFINE_TRANSPORT": GEN / "p2415_s1365_phase_frequency_affine_transport_nonautomorphism_certificate.json",
    "SCRATCH_PHASE_ZERO_CELL_SIGN": ROOT / "scratch" / "bridge_strict_completion_phase_zero_cell_sign_certificate_report.md",
}

STRICT_OMEGA = Fraction(743, 4000)
STRICT_PHI = Fraction(13, 80)
DOMAIN = list(range(12))
GRID_FACTORS = [Fraction(-1, 1), Fraction(-1, 2), Fraction(0, 1), Fraction(1, 2), Fraction(1, 1)]
NEGATIVE_EXPORT_FLAGS = [
    "finite_sign_pattern_phase_frequency_source_exported",
    "strict_phase_frequency_source_exported",
    "strict_dynamical_source_for_A_P_D_exported",
    "strict_damping_beta_eta_source_exported",
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
        "new_packet": "P2564|S1514|phase frequency finite sign cell|finite sign-cell nonidentifiability|phase sign cell nonidentifiability",
        "intended_research_nonduplication": "sign-cell|sign cell nonidentifiability|phase.*open cell|finite sign.*phase|sign pattern.*omega.*phi|omega.*phi.*sign pattern",
        "precursor_phase_signs": "phase zero cell sign|finite sign pattern|phase-sign GF2|phase sign z2|phase/frequency affine transport",
        "immediate_precursors": "P2563|S1513|rational winding quotient|P2561|S1511|strict_phase_frequency_source|P2415|S1365",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def file_presence(path: Any) -> dict[str, Any]:
    if hasattr(path, "exists") and path.exists():
        text = path.read_text(encoding="utf-8")
        return {"path": rel(path), "present": True, "sha256": hashlib.sha256(text.encode()).hexdigest(), "bytes": len(text.encode())}
    return {"path": rel(path), "present": False}


def sign(value: float) -> int:
    if value > 0.0:
        return 1
    if value < 0.0:
        return -1
    return 0


def zero_clearance(theta: float) -> float:
    nearest = min(abs(theta - (math.pi / 2.0 + k * math.pi)) for k in range(-8, 9))
    return nearest


def strict_phase_rows() -> list[dict[str, Any]]:
    rows = []
    for d in DOMAIN:
        theta = float(STRICT_OMEGA) * d + float(STRICT_PHI)
        rows.append({
            "d": d,
            "theta": theta,
            "cos_theta": math.cos(theta),
            "sign": sign(math.cos(theta)),
            "zero_clearance_radians": zero_clearance(theta),
        })
    return rows


def perturbation_box(rows: list[dict[str, Any]]) -> dict[str, Any]:
    min_row = min(rows, key=lambda row: row["zero_clearance_radians"])
    max_d = max(DOMAIN)
    epsilon = min_row["zero_clearance_radians"] / (2.0 * (max_d + 1.0))
    return {
        "min_clearance_row": min_row,
        "max_d": max_d,
        "epsilon_omega": epsilon,
        "epsilon_phi": epsilon,
        "box_area": (2.0 * epsilon) ** 2,
        "proof_bound": "if |delta_omega|<=epsilon and |delta_phi|<=epsilon then |d*delta_omega+delta_phi| <= (max_d+1)*epsilon = min_clearance/2, so no audited node crosses a cosine zero",
    }


def signs_for(omega: float, phi: float) -> list[int]:
    return [sign(math.cos(omega * d + phi)) for d in DOMAIN]


def grid_witnesses(box: dict[str, Any], strict_signs: list[int]) -> list[dict[str, Any]]:
    eps_o = box["epsilon_omega"]
    eps_p = box["epsilon_phi"]
    rows = []
    for omega_factor in GRID_FACTORS:
        for phi_factor in GRID_FACTORS:
            delta_omega = float(omega_factor) * eps_o
            delta_phi = float(phi_factor) * eps_p
            omega = float(STRICT_OMEGA) + delta_omega
            phi = float(STRICT_PHI) + delta_phi
            signs = signs_for(omega, phi)
            rows.append({
                "omega_factor": str(omega_factor),
                "phi_factor": str(phi_factor),
                "delta_omega": delta_omega,
                "delta_phi": delta_phi,
                "omega": omega,
                "phi": phi,
                "same_sign_pattern_as_strict": signs == strict_signs,
                "tuple_differs_from_strict": delta_omega != 0.0 or delta_phi != 0.0,
                "signs": signs,
            })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2563_payload = load_json(SOURCE_FILES["P2563_RATIONAL_WINDING_OBSTRUCTION"])
    p2561_payload = load_json(SOURCE_FILES["P2561_POST_DAMPING_RESIDUAL_BRIDGE"])
    p2415_payload = load_json(SOURCE_FILES["P2415_PHASE_AFFINE_TRANSPORT"])
    p2563 = theorem(p2563_payload, "phase_frequency_rational_winding_quotient_obstruction_certificate")
    p2561 = theorem(p2561_payload, "strict_completion_post_damping_residual_bridge_two_key_certificate")
    p2415 = theorem(p2415_payload, "phase_frequency_affine_transport_nonautomorphism_certificate")
    phase_rows = strict_phase_rows()
    strict_signs = [row["sign"] for row in phase_rows]
    box = perturbation_box(phase_rows)
    witnesses = grid_witnesses(box, strict_signs)
    nontrivial_witnesses = [row for row in witnesses if row["tuple_differs_from_strict"] and row["same_sign_pattern_as_strict"]]
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2564_T1_phase_frequency_finite_sign_cell_nonidentifiability_certificate",
        "audited_chain": ["P2563/S1513", "P2561/S1511", "P2415/S1365", "scratch phase-zero cell-sign"],
        "frontier_atom_under_attack": "strict_phase_frequency_source",
        "source_class_audited": "finite d=0..11 phase-sign/topological-bit pattern as selector for exact omega and phi",
        "p2563_rational_winding_obstruction_inherited": p2563.get("exact_rational_winding_quotient_source_rejected") is True,
        "p2561_phase_frequency_residual_atom_inherited": "strict_phase_frequency_source" in p2561.get("residual_atoms_after_hypothetical_damping_source", []),
        "p2415_affine_nonautomorphism_inherited": p2415.get("no_z12_unit_offset_reindex_matches_strict_sign_pattern") is True,
        "strict_phase_rows_d0_to_d11": phase_rows,
        "strict_sign_pattern_d0_to_d11": strict_signs,
        "sign_pattern_positive_count": sum(1 for value in strict_signs if value > 0),
        "sign_pattern_negative_count": sum(1 for value in strict_signs if value < 0),
        "certified_open_sign_cell_box": box,
        "grid_factor_set": [str(value) for value in GRID_FACTORS],
        "grid_witness_count": len(witnesses),
        "grid_witnesses_preserving_sign_count": sum(1 for row in witnesses if row["same_sign_pattern_as_strict"]),
        "nontrivial_same_sign_witness_count": len(nontrivial_witnesses),
        "sample_nontrivial_same_sign_witnesses": nontrivial_witnesses[:8],
        "finite_sign_pattern_selects_unique_phase_frequency_tuple": False,
        "finite_sign_pattern_has_open_continuum_of_phase_frequency_realizations": box["box_area"] > 0.0 and len(nontrivial_witnesses) > 0,
        "required_new_source_to_escape_obstruction": "a metric/dynamical/selector principle stronger than the finite phase-sign/topological-bit pattern, capable of selecting omega=743/4000 and phi=13/80 inside the certified open sign cell",
        "recommended_next_honest_step": (
            "Do not treat GF(2) or finite phase-sign reconstruction as a source for the exact strict phase/frequency tuple. The next honest step is to add a stronger selector candidate with a computable objective inside the open sign cell, then test whether it has a unique minimizer at omega=743/4000 and phi=13/80 without discharging QW-2191 by assumption."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2563_inherited": theorem_export["p2563_rational_winding_obstruction_inherited"],
        "p2561_phase_frequency_atom_inherited": theorem_export["p2561_phase_frequency_residual_atom_inherited"],
        "p2415_nonautomorphism_inherited": theorem_export["p2415_affine_nonautomorphism_inherited"],
        "finite_phase_rows_are_12": len(theorem_export["strict_phase_rows_d0_to_d11"]) == 12,
        "open_box_has_positive_area": theorem_export["certified_open_sign_cell_box"]["box_area"] > 0.0,
        "all_grid_witnesses_preserve_signs": theorem_export["grid_witnesses_preserving_sign_count"] == theorem_export["grid_witness_count"],
        "nontrivial_same_sign_witnesses_exist": theorem_export["nontrivial_same_sign_witness_count"] > 0,
        "finite_sign_pattern_not_unique": theorem_export["finite_sign_pattern_selects_unique_phase_frequency_tuple"] is False,
        "phase_frequency_source_not_exported": theorem_export["strict_phase_frequency_source_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2564",
        "stage_id": "S1514",
        "status": "P2564_PHASE_FREQUENCY_FINITE_SIGN_CELL_NONIDENTIFIABILITY_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "phase_frequency_finite_sign_cell_nonidentifiability_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2563_RATIONAL_WINDING_OBSTRUCTION": sha256_json(p2563_payload),
                "P2561_POST_DAMPING_RESIDUAL_BRIDGE": sha256_json(p2561_payload),
                "P2415_PHASE_AFFINE_TRANSPORT": sha256_json(p2415_payload),
                "SCRATCH_PHASE_ZERO_CELL_SIGN": file_presence(SOURCE_FILES["SCRATCH_PHASE_ZERO_CELL_SIGN"]),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["phase_frequency_finite_sign_cell_nonidentifiability_certificate"]["theorem_export"]
    box = t["certified_open_sign_cell_box"]
    lines = [
        "# P2564/S1514 phase/frequency finite sign-cell nonidentifiability certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Source class audited: `{t['source_class_audited']}`.",
        f"- P2563 rational-winding obstruction inherited: `{t['p2563_rational_winding_obstruction_inherited']}`.",
        f"- P2561 phase/frequency residual atom inherited: `{t['p2561_phase_frequency_residual_atom_inherited']}`.",
        f"- P2415 affine nonautomorphism inherited: `{t['p2415_affine_nonautomorphism_inherited']}`.",
        f"- Strict sign pattern d=0..11: `{t['strict_sign_pattern_d0_to_d11']}`.",
        f"- Certified epsilon omega: `{box['epsilon_omega']}`.",
        f"- Certified epsilon phi: `{box['epsilon_phi']}`.",
        f"- Open sign-cell box area: `{box['box_area']}`.",
        f"- Grid witnesses preserving signs: `{t['grid_witnesses_preserving_sign_count']}/{t['grid_witness_count']}`.",
        f"- Nontrivial same-sign witnesses: `{t['nontrivial_same_sign_witness_count']}`.",
        f"- Finite sign pattern selects unique tuple: `{t['finite_sign_pattern_selects_unique_phase_frequency_tuple']}`.", "",
        "## Proof interpretation", "",
        "The strict finite sign pattern has a positive clearance from every audited cosine zero.  A small open box around `(omega, phi)=(743/4000,13/80)` therefore preserves all audited signs, so finite GF(2)/phase-sign data cannot by itself select the exact strict phase/frequency tuple.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No finite-sign phase/frequency source, strict phase/frequency source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['phase_frequency_finite_sign_cell_nonidentifiability_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2564/S1514` audits whether the finite strict phase-sign/topological-bit pattern on `d=0..11` can itself source the exact strict phase/frequency tuple.  It cannot: the strict tuple has positive clearance from all audited cosine-zero boundaries, yielding an explicit open box around `omega=743/4000`, `phi=13/80` that preserves the same finite sign pattern.  A 25-point perturbation grid inside the box preserves all signs, so finite sign/GF(2) reconstruction is not a unique phase/frequency selector.
""".strip()
    lag_section = """
`P2564/S1514` blocks promotion of finite phase-sign data into role-bearing phase/frequency dynamics in `L_total`.  The exact tuple still requires a stronger strict selector/source principle inside the open sign cell; the certificate exports no bridge theorem, role-transfer theorem, QW-2191 discharge, or ToE closure.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2564/S1514 phase/frequency finite sign-cell nonidentifiability guard", "## P2564/S1514 phase/frequency finite sign-cell nonidentifiability guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2564/S1514 phase/frequency finite sign-cell nonidentifiability Ltotal guard", "## P2564/S1514 phase/frequency finite sign-cell nonidentifiability Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
