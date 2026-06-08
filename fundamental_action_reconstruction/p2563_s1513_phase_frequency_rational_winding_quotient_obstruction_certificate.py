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
OUT = GEN / "p2563_s1513_phase_frequency_rational_winding_quotient_obstruction_certificate.json"
MD = GEN / "p2563_s1513_phase_frequency_rational_winding_quotient_obstruction_certificate.md"

SOURCE_FILES = {
    "P2562_SHORTCUT_NONPROMOTION": GEN / "p2562_s1512_physical_ontology_shortcut_nonpromotion_audit.json",
    "P2561_POST_DAMPING_RESIDUAL_BRIDGE": GEN / "p2561_s1511_strict_completion_post_damping_residual_bridge_two_key_certificate.json",
    "P2415_PHASE_AFFINE_TRANSPORT": GEN / "p2415_s1365_phase_frequency_affine_transport_nonautomorphism_certificate.json",
    "F326_PHASE_FREQUENCY_OBSTRUCTION": ROOT / "F326_FIRST_ACTUAL_LEGACY_TO_STRICT_KERNEL_PHASE_FREQUENCY_NONCONFORMAL_OBSTRUCTION_WITNESS_PACKET.md",
}

STRICT_OMEGA = Fraction(743, 4000)
STRICT_PHI = Fraction(13, 80)
LEGACY_OMEGA_PI_COEFFICIENT = Fraction(1, 4)
LEGACY_PHI_PI_COEFFICIENT = Fraction(1, 6)
SEARCH_NUMERATOR_BOUND = 96
SEARCH_DENOMINATOR_BOUND = 96
DOMAIN = list(range(12))
NEGATIVE_EXPORT_FLAGS = [
    "rational_winding_quotient_phase_frequency_source_exported",
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
        "new_packet": "P2563|S1513|rational winding quotient|winding quotient obstruction|phase-frequency rational winding|phase/frequency rational-winding",
        "intended_research_nonduplication": "rational.*winding.*phase|phase.*rational.*winding|pi-rational phase|rational pi multiple|Diophantine.*phase/frequency|phase/frequency.*Diophantine",
        "precursor_phase_frequency": "P2415|S1365|phase frequency affine transport|phase-frequency affine transport nonautomorphism|F326|phase/frequency nonconformal obstruction",
        "immediate_precursors": "P2562|S1512|physical ontology shortcut|P2561|S1511|post-damping residual bridge|strict_phase_frequency_source",
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


def rational_pi_multiple_search(target: Fraction) -> dict[str, Any]:
    target_float = float(target)
    best: dict[str, Any] | None = None
    for denominator in range(1, SEARCH_DENOMINATOR_BOUND + 1):
        for numerator in range(-SEARCH_NUMERATOR_BOUND, SEARCH_NUMERATOR_BOUND + 1):
            coefficient = Fraction(numerator, denominator)
            value = float(coefficient) * math.pi
            residual = abs(value - target_float)
            row = {
                "pi_coefficient": str(coefficient),
                "numerator": numerator,
                "denominator": denominator,
                "candidate_float": value,
                "target_float": target_float,
                "abs_residual": residual,
            }
            if best is None or residual < best["abs_residual"]:
                best = row
    assert best is not None
    return best


def exact_obstruction_rows() -> list[dict[str, Any]]:
    return [
        {
            "slot": "omega",
            "legacy_pi_coefficient": str(LEGACY_OMEGA_PI_COEFFICIENT),
            "strict_value": str(STRICT_OMEGA),
            "allowed_winding_quotient_form": "q*pi for rational q",
            "exact_equality_possible_without_new_pi_cancelling_source": False,
            "proof_reason": "Every nonzero rational multiple of pi is irrational, while 743/4000 is rational and nonzero; q=0 gives 0, not 743/4000.",
            "bounded_best_pi_multiple": rational_pi_multiple_search(STRICT_OMEGA),
        },
        {
            "slot": "phi",
            "legacy_pi_coefficient": str(LEGACY_PHI_PI_COEFFICIENT),
            "strict_value": str(STRICT_PHI),
            "allowed_winding_quotient_form": "q*pi for rational q",
            "exact_equality_possible_without_new_pi_cancelling_source": False,
            "proof_reason": "Every nonzero rational multiple of pi is irrational, while 13/80 is rational and nonzero; q=0 gives 0, not 13/80.",
            "bounded_best_pi_multiple": rational_pi_multiple_search(STRICT_PHI),
        },
    ]


def finite_phase_table() -> list[dict[str, Any]]:
    rows = []
    for d in DOMAIN:
        legacy_phase = math.pi * (float(LEGACY_OMEGA_PI_COEFFICIENT) * d + float(LEGACY_PHI_PI_COEFFICIENT))
        strict_phase = float(STRICT_OMEGA) * d + float(STRICT_PHI)
        rows.append({
            "d": d,
            "legacy_phase": legacy_phase,
            "strict_phase": strict_phase,
            "phase_gap_abs": abs(legacy_phase - strict_phase),
            "legacy_cos": math.cos(legacy_phase),
            "strict_cos": math.cos(strict_phase),
            "cos_gap_abs": abs(math.cos(legacy_phase) - math.cos(strict_phase)),
            "legacy_sign": 1 if math.cos(legacy_phase) > 0 else -1 if math.cos(legacy_phase) < 0 else 0,
            "strict_sign": 1 if math.cos(strict_phase) > 0 else -1 if math.cos(strict_phase) < 0 else 0,
        })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2562_payload = load_json(SOURCE_FILES["P2562_SHORTCUT_NONPROMOTION"])
    p2561_payload = load_json(SOURCE_FILES["P2561_POST_DAMPING_RESIDUAL_BRIDGE"])
    p2415_payload = load_json(SOURCE_FILES["P2415_PHASE_AFFINE_TRANSPORT"])
    p2562 = theorem(p2562_payload, "physical_ontology_shortcut_nonpromotion_audit")
    p2561 = theorem(p2561_payload, "strict_completion_post_damping_residual_bridge_two_key_certificate")
    p2415 = theorem(p2415_payload, "phase_frequency_affine_transport_nonautomorphism_certificate")
    obstruction = exact_obstruction_rows()
    phase_rows = finite_phase_table()
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2563_T1_phase_frequency_rational_winding_quotient_obstruction_certificate",
        "audited_chain": ["P2562/S1512", "P2561/S1511", "P2415/S1365", "F326"],
        "frontier_atom_under_attack": "strict_phase_frequency_source",
        "p2562_shortcut_nonpromotion_inherited": p2562.get("physical_ontology_document_can_replace_source_theorem") is False,
        "p2561_phase_frequency_residual_atom_inherited": "strict_phase_frequency_source" in p2561.get("residual_atoms_after_hypothetical_damping_source", []),
        "p2415_affine_nonautomorphism_inherited": p2415.get("no_z12_unit_offset_reindex_matches_strict_sign_pattern") is True,
        "allowed_source_class_audited": "pure rational winding/cycle quotient of legacy pi-rational phase data",
        "exact_obstruction_rows": obstruction,
        "exact_rational_winding_quotient_source_rejected": all(not row["exact_equality_possible_without_new_pi_cancelling_source"] for row in obstruction),
        "bounded_pi_multiple_search": {
            "numerator_bound": SEARCH_NUMERATOR_BOUND,
            "denominator_bound": SEARCH_DENOMINATOR_BOUND,
            "omega_best": obstruction[0]["bounded_best_pi_multiple"],
            "phi_best": obstruction[1]["bounded_best_pi_multiple"],
        },
        "finite_phase_table_d0_to_d11": phase_rows,
        "finite_domain_size": len(phase_rows),
        "max_phase_gap_abs_d0_to_d11": max(row["phase_gap_abs"] for row in phase_rows),
        "max_cos_gap_abs_d0_to_d11": max(row["cos_gap_abs"] for row in phase_rows),
        "sign_mismatch_count_d0_to_d11": sum(1 for row in phase_rows if row["legacy_sign"] != row["strict_sign"]),
        "required_new_source_to_escape_obstruction": "a non-winding, non-pure-rational-quotient strict source that introduces the rational strict phase/frequency tuple or a justified pi-cancelling renormalization map",
        "recommended_next_honest_step": (
            "Do not try another pure winding/cycle quotient for phase/frequency. The next honest step is to formulate an explicit non-winding strict phase/frequency source theorem: a nadsoliton dynamical equation or selector premise that exports omega=743/4000 and phi=13/80 directly, while keeping QW-2191 open unless that source also breaks the selector symmetry."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2562_inherited": theorem_export["p2562_shortcut_nonpromotion_inherited"],
        "p2561_phase_frequency_atom_inherited": theorem_export["p2561_phase_frequency_residual_atom_inherited"],
        "p2415_nonautomorphism_inherited": theorem_export["p2415_affine_nonautomorphism_inherited"],
        "two_exact_obstruction_rows": len(theorem_export["exact_obstruction_rows"]) == 2,
        "rational_winding_source_rejected": theorem_export["exact_rational_winding_quotient_source_rejected"],
        "bounded_search_bounds_recorded": theorem_export["bounded_pi_multiple_search"]["numerator_bound"] == SEARCH_NUMERATOR_BOUND and theorem_export["bounded_pi_multiple_search"]["denominator_bound"] == SEARCH_DENOMINATOR_BOUND,
        "finite_phase_table_has_12_rows": theorem_export["finite_domain_size"] == 12,
        "phase_frequency_source_not_exported": theorem_export["strict_phase_frequency_source_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2563",
        "stage_id": "S1513",
        "status": "P2563_PHASE_FREQUENCY_RATIONAL_WINDING_QUOTIENT_OBSTRUCTION_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "phase_frequency_rational_winding_quotient_obstruction_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2562_SHORTCUT_NONPROMOTION": sha256_json(p2562_payload),
                "P2561_POST_DAMPING_RESIDUAL_BRIDGE": sha256_json(p2561_payload),
                "P2415_PHASE_AFFINE_TRANSPORT": sha256_json(p2415_payload),
                "F326_PHASE_FREQUENCY_OBSTRUCTION": file_presence(SOURCE_FILES["F326_PHASE_FREQUENCY_OBSTRUCTION"]),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["phase_frequency_rational_winding_quotient_obstruction_certificate"]["theorem_export"]
    search = t["bounded_pi_multiple_search"]
    lines = [
        "# P2563/S1513 phase/frequency rational-winding quotient obstruction certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2562 shortcut nonpromotion inherited: `{t['p2562_shortcut_nonpromotion_inherited']}`.",
        f"- P2561 phase/frequency residual atom inherited: `{t['p2561_phase_frequency_residual_atom_inherited']}`.",
        f"- P2415 affine nonautomorphism inherited: `{t['p2415_affine_nonautomorphism_inherited']}`.",
        f"- Allowed source class audited: `{t['allowed_source_class_audited']}`.",
        f"- Exact rational-winding quotient source rejected: `{t['exact_rational_winding_quotient_source_rejected']}`.",
        f"- Bounded pi-multiple search bounds: `|numerator| <= {search['numerator_bound']}`, `denominator <= {search['denominator_bound']}`.",
        f"- Best omega residual in bounded search: `{search['omega_best']['abs_residual']}` at coefficient `{search['omega_best']['pi_coefficient']}`.",
        f"- Best phi residual in bounded search: `{search['phi_best']['abs_residual']}` at coefficient `{search['phi_best']['pi_coefficient']}`.",
        f"- Finite phase rows d=0..11: `{t['finite_domain_size']}`.",
        f"- Sign mismatches d=0..11: `{t['sign_mismatch_count_d0_to_d11']}`.", "",
        "## Proof interpretation", "",
        "A pure rational winding/cycle quotient of the legacy phase layer stays in the class `q*pi` with rational `q`. The strict targets `743/4000` and `13/80` are nonzero rationals. Therefore exact equality would require a new pi-cancelling or non-winding strict source, not just topological quotient bookkeeping.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No rational-winding phase/frequency source, strict phase/frequency source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['phase_frequency_rational_winding_quotient_obstruction_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2563/S1513` attacks the phase/frequency residual atom from P2561 by auditing a precise source class: pure rational winding/cycle quotients of the legacy `pi/4, pi/6` phase data.  The proof obstruction is exact: such quotients remain rational multiples of `pi`, while the strict targets `omega=743/4000` and `phi=13/80` are nonzero rationals.  A bounded search over `|numerator| <= 96`, `denominator <= 96` records best numerical approximants, but does not remove the exact obstruction.  Therefore phase/frequency still needs a non-winding strict source or an explicitly justified pi-cancelling renormalization map.
""".strip()
    lag_section = """
`P2563/S1513` blocks promotion of a pure topological winding quotient into role-bearing phase/frequency dynamics in `L_total`.  Phase/frequency remains a strict source obligation; the certificate exports no bridge theorem, role-transfer theorem, QW-2191 discharge, or ToE closure.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2563/S1513 phase/frequency rational-winding quotient obstruction guard", "## P2563/S1513 phase/frequency rational-winding quotient obstruction guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2563/S1513 phase/frequency rational-winding quotient obstruction Ltotal guard", "## P2563/S1513 phase/frequency rational-winding quotient obstruction Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
