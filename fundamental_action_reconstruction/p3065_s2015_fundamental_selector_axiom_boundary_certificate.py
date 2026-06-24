#!/usr/bin/env python3
"""P3065/S2015: fundamental selector axiom boundary certificate.

The user asks to stop letting the non-premise selector search block the study of
an informational nadsoliton theory.  P3065 therefore constructs the missing
axiom-boundary object honestly: a two-member family of axiom-augmented theories
T_sigma, "the universe with a given orientation choice".  This does not derive
sigma_selector, does not discharge QW-2191, and does not become strict selector
closure.  It gives a finite, computable contract for using sigma_selector as a
boundary parameter while keeping strict and axiom-augmented exports separated.
"""
from __future__ import annotations

import hashlib, json, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3064_s2014_strict_polarity_source_theorem_sat import OUT as P3064

OUT = GEN / "p3065_s2015_fundamental_selector_axiom_boundary_certificate.json"
MD = GEN / "p3065_s2015_fundamental_selector_axiom_boundary_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "selector_as_axiom_boundary_parameter": r"selector.*axiom|selection.*axiom|sigma_selector.*parameter|boundary parameter.*selector|given orientation",
    "axiom_augmented_not_strict_selector_closure": r"axiom[-_ ]augmented.*selector|non-strict.*selector|QW-2191.*axiom-free|explicit.*selection.*axiom",
    "universe_with_given_orientation_family": r"orientation.*choice|given.*orientation|theory.*orientation|T_sigma|sigma.*orientation",
}

SIGMA_CHOICES = ("sigma_plus", "sigma_minus")
STRICT_PRIMITIVE_ATOMS = ("strict_source_plus", "strict_source_minus", "strict_coupling_identity", "strict_coupling_flip")
CURRENT_STRICT_EXPORTED_ATOMS: tuple[str, ...] = ()


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:10]})
    return rows


def evaluate_theory_row(sigma_axiom: str | None, strict_atoms: tuple[str, ...]) -> dict[str, Any]:
    strict_set = set(strict_atoms)
    source_count = len(strict_set.intersection({"strict_source_plus", "strict_source_minus"}))
    coupling_count = len(strict_set.intersection({"strict_coupling_identity", "strict_coupling_flip"}))
    strict_selector_export = source_count == 1 and coupling_count == 1
    axiom_augmented_admitted = sigma_axiom in SIGMA_CHOICES
    branch_orientation = {"sigma_plus": "+omega/G_plus", "sigma_minus": "-omega/G_minus"}.get(sigma_axiom)
    blocked_exports = []
    if axiom_augmented_admitted and not strict_selector_export:
        blocked_exports.extend(["strict_selector_closure", "QW-2191_discharge", "derived_sigma_selector", "role_bearing_L_total", "ToE_closure"])
    if not axiom_augmented_admitted:
        blocked_exports.append("orientation_conditioned_branch_use")
    return {
        "sigma_axiom": sigma_axiom,
        "strict_atoms": list(strict_atoms),
        "axiom_augmented_theory_admitted": axiom_augmented_admitted,
        "branch_orientation": branch_orientation,
        "strict_selector_export": strict_selector_export,
        "sigma_status": "fundamental_boundary_parameter_not_derivation" if axiom_augmented_admitted else "absent",
        "may_continue_nadsoliton_analysis_conditioned_on_sigma": axiom_augmented_admitted,
        "blocked_exports": blocked_exports,
    }


def build_payload() -> dict[str, Any]:
    p3064 = read_json(P3064)
    grep_rows = content_grep()
    rows = [evaluate_theory_row(sigma, CURRENT_STRICT_EXPORTED_ATOMS) for sigma in (None,) + SIGMA_CHOICES]
    axiom_rows = [r for r in rows if r["axiom_augmented_theory_admitted"]]
    proof_obligations = [
        {"obligation": "content_first_grep_before_axiom_boundary", "satisfied": True, "detail": "searched by selector axiom/boundary-parameter/orientation-choice content, not by P-numbers only"},
        {"obligation": "construct_fundamental_selector_axiom_object", "satisfied": True, "detail": "constructed T_sigma as the two-row family sigma_plus/sigma_minus"},
        {"obligation": "separate_axiom_augmented_use_from_strict_derivation", "satisfied": True, "detail": "axiom rows admit conditioned branch work but keep strict_selector_export false"},
        {"obligation": "resume_information_nadsoliton_work_without_false_selector_claim", "satisfied": True, "detail": "orientation-conditioned work is allowed with sigma recorded as boundary datum"},
        {"obligation": "derive_sigma_selector_nonpremise", "satisfied": False, "detail": "P3065 explicitly does not derive the selector and does not fix the P3064 primitive atoms"},
    ]
    return {
        "status": "P3065_FUNDAMENTAL_SELECTOR_AXIOM_BOUNDARY_CERTIFICATE_AXIOM_AUGMENTED_ONLY",
        "input_hashes": {"P3064": hashlib.sha256(P3064.read_bytes()).hexdigest() if P3064.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "fundamental_selector_axiom": {
                "object": "FundamentalSelectorAxiomBoundary",
                "formal_family": "T_sigma = strict_nadsoliton_information_theory + A_selector(sigma), sigma in {+1,-1}",
                "interpretation": "theory = universe with a given orientation choice; sigma_selector is a boundary parameter, not a derivation",
                "allowed_use": "condition all downstream orientation-sensitive constructions on the recorded sigma value",
                "forbidden_use": "do not rebrand the axiom as non-premise strict selector source, observed physics, L_total, bridge closure, role transfer, or ToE",
            },
            "axiom_boundary_rows": rows,
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(r["hit_count"] for r in grep_rows),
            "sigma_choices": len(SIGMA_CHOICES),
            "boundary_rows": len(rows),
            "axiom_augmented_admitted_rows": len(axiom_rows),
            "strict_selector_export_rows": sum(1 for r in rows if r["strict_selector_export"]),
            "current_strict_exported_atoms": len(CURRENT_STRICT_EXPORTED_ATOMS),
            "conditioned_nadsoliton_continuation_rows": sum(1 for r in rows if r["may_continue_nadsoliton_analysis_conditioned_on_sigma"]),
            "p3064_current_row_accepted": p3064.get("finite_certificate", {}).get("current_row_accepted"),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(1 for o in proof_obligations if o["satisfied"]),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "axiom_boundary_result": "P3065 accepts sigma_selector only as an explicit fundamental boundary axiom.  The finite boundary table has two admitted axiom-augmented theories, T_sigma_plus and T_sigma_minus, and zero strict selector-export rows.  This unblocks orientation-conditioned informational-nadsoliton analysis while preserving the P3064 no-current-strict-export result.",
            "negative_export_flags": {k: False for k in ["sigma_selector_derived", "p3064_primitive_atom_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"axiom_augmented_T_sigma_family_constructed": True, "sigma_selector_boundary_parameter_recorded": True, "orientation_conditioned_nadsoliton_analysis_admitted": True},
            "next_honest_step": "Proceed under the explicit T_sigma axiom boundary: choose one orientation branch as an assumed boundary datum for the next computation, then build a finite orientation-conditioned informational-nadsoliton object that does not claim selector derivation.  The best next proof-grade move is a sigma-conditioned state-transition/observable-invariant matrix showing which downstream constructions depend on sigma and which are sigma-invariant; keep all exports marked axiom-augmented unless a future strict source theorem supplies the missing P3064 atoms.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3065/S2015 fundamental selector axiom boundary certificate", "", f"Status: `{payload['status']}`", "",
        "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- sigma choices: `{c['sigma_choices']}`",
        f"- boundary rows: `{c['boundary_rows']}`",
        f"- axiom-augmented admitted rows: `{c['axiom_augmented_admitted_rows']}`",
        f"- strict selector-export rows: `{c['strict_selector_export_rows']}`",
        f"- conditioned nadsoliton continuation rows: `{c['conditioned_nadsoliton_continuation_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        "", "## Decision", payload["decision"]["axiom_boundary_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3065/S2015 fundamental selector axiom boundary certificate", "## P3065/S2015 fundamental selector axiom boundary certificate\n\n`P3065/S2015` constructs the explicit axiom-boundary object requested after `P3064`: `T_sigma = strict_nadsoliton_information_theory + A_selector(sigma)`, with `sigma in {+1,-1}`.  The finite boundary table has `3` rows (`sigma` absent, `sigma_plus`, `sigma_minus`), `2` axiom-augmented admitted branches, and `0` strict selector-export rows.  Thus `sigma_selector` may be used as a recorded boundary parameter for orientation-conditioned informational-nadsoliton analysis, but it remains non-derived and does not discharge `QW-2191`, export strict selector closure, create observed physics, `L_total`, bridge/role transfer, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3065/S2015 selector axiom is not an `L_total` source", "## P3065/S2015 selector axiom is not an `L_total` source\n\n`P3065/S2015` adds no variational term.  It records an axiom-augmented orientation boundary condition for branch-conditioned study; no unit-bearing action density, EOM term, Hamiltonian term, or nonproxy coupling is exported.\n")
    append_once(AGENTS, "Current fundamental selector axiom boundary guardrail (P3065/S2015, 2026-06-24)", "## Current fundamental selector axiom boundary guardrail (P3065/S2015, 2026-06-24)\n\n- P3065 accepts `sigma_selector` only as an explicit fundamental boundary axiom: `T_sigma = strict_nadsoliton_information_theory + A_selector(sigma)`, with `sigma in {+1,-1}`; this means theory = universe with a given orientation choice.\n- The finite boundary table has `2` admitted axiom-augmented branches and `0` strict selector-export rows.  `sigma_selector` is therefore a recorded boundary parameter, not a non-premise derivation or strict source theorem.\n- This axiom boundary may unblock orientation-conditioned study of the informational nadsoliton, but every dependent result must remain marked axiom-augmented/conditioned-on-`sigma` unless a later strict theorem supplies the missing P3064 source/coupling atoms.\n- Do not promote the axiom boundary to observed physics, `QW-2191` discharge, strict selector closure, `L_total`, bridge/role-transfer, or ToE closure.  The next proof-grade move should build a finite `sigma`-conditioned state-transition/observable-invariant matrix separating sigma-dependent from sigma-invariant downstream objects.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
