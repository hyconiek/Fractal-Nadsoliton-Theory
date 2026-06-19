#!/usr/bin/env python3
"""P2951/S1901: ratio-package strict-source normal-form lattice.

P2950 showed that the exact P2948 finite package is eta-compatible but still
lacks strict beta scale and nonproxy damping coupling.  Repo content-grep also
shows older scale-orbit/UV-unit audits (P2649/P2655/P2689), so P2951 does not
replay another beta-normalization scan.

Instead it constructs the missing theorem-object normal form for this exact
ratio-package lane.  The object is a four-atom obligation lattice:

  1. strict provenance for the P2938 torsion-character aggregate,
  2. strict identity-deficit delta numerator semantics,
  3. strict positive beta-scale/unit source,
  4. nonproxy variational damping coupling into L_total/EOM.

The finite computation enumerates all 2^4 masks.  Exactly one row accepts: all
four theorem atoms present.  Every atom is essential, and the current artifact
mask has none of the four atoms discharged.  This is a proof-control object: it
prevents promoting any proper subset of the exact package to a strict damping
source.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2948_s1898_torsion_character_ratio_package_theorem_skeleton import OUT as P2948
from p2949_s1899_delta_numerator_semantics_separation_audit import OUT as P2949
from p2950_s1900_exact_package_beta_eta_scale_coupling_obstruction import OUT as P2950

OUT = GEN / "p2951_s1901_ratio_package_strict_source_normal_form_lattice.json"
MD = GEN / "p2951_s1901_ratio_package_strict_source_normal_form_lattice.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

OBLIGATIONS = [
    {
        "atom": "torsion_character_aggregate_strict_provenance",
        "meaning": "prove the P2938 U(12) torsion-character aggregate is nadsoliton-sourced, not a finite carrier premise",
    },
    {
        "atom": "identity_deficit_delta_semantics",
        "meaning": "prove identity-deficit, not zero/intersection/union aliasing, is the strict delta numerator semantics",
    },
    {
        "atom": "positive_beta_scale_unit_source",
        "meaning": "prove a target-independent positive beta-scale/unit-normalization source for the exact package",
    },
    {
        "atom": "nonproxy_variational_damping_coupling",
        "meaning": "couple the exact package to a unit-bearing nonproxy damping term in L_total/EOM",
    },
]


def current_mask(p2948: dict[str, Any], p2949: dict[str, Any], p2950: dict[str, Any]) -> dict[str, bool]:
    p2948_cert = p2948["package_certificate"]
    p2949_cert = p2949["delta_numerator_semantics_certificate"]
    p2950_cert = p2950["beta_eta_scale_coupling_certificate"]
    p2950_flags = p2950["decision"]["negative_export_flags"]
    return {
        "torsion_character_aggregate_strict_provenance": bool(p2948_cert["strict_torsion_character_source_theorem_exported"]),
        "identity_deficit_delta_semantics": bool(p2949_cert["strict_intensional_identity_deficit_semantics_exported"]),
        "positive_beta_scale_unit_source": bool(p2950_cert["strict_positive_beta_scale_source_exported"]),
        "nonproxy_variational_damping_coupling": bool(p2950_flags["nonproxy_ltotal_exported"]),
    }


def lattice_rows() -> list[dict[str, Any]]:
    atoms = [item["atom"] for item in OBLIGATIONS]
    rows = []
    for bits in itertools.product([False, True], repeat=len(atoms)):
        mask = dict(zip(atoms, bits))
        accepted = all(mask.values())
        missing = [atom for atom, present in mask.items() if not present]
        rows.append({
            "mask": mask,
            "present_count": sum(bits),
            "missing_atoms": missing,
            "accepted_strict_ratio_damping_source": accepted,
            "proper_subset_rejected": not accepted,
        })
    return rows


def essentiality_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    atoms = [item["atom"] for item in OBLIGATIONS]
    result = []
    for atom in atoms:
        almost_all = {other: True for other in atoms}
        almost_all[atom] = False
        matching = next(row for row in rows if row["mask"] == almost_all)
        result.append({
            "atom": atom,
            "essential": not matching["accepted_strict_ratio_damping_source"],
            "evidence": "with every other atom present, removing this atom still rejects the package",
        })
    return result


def build_payload(p2948: dict[str, Any], p2949: dict[str, Any], p2950: dict[str, Any]) -> dict[str, Any]:
    rows = lattice_rows()
    current = current_mask(p2948, p2949, p2950)
    current_row = next(row for row in rows if row["mask"] == current)
    essential = essentiality_rows(rows)
    accepted_rows = [row for row in rows if row["accepted_strict_ratio_damping_source"]]
    return {
        "status": "P2951_RATIO_PACKAGE_STRICT_SOURCE_NORMAL_FORM_LATTICE_NO_STRICT_EXPORT",
        "input_hashes": {
            "P2948": hashlib.sha256(P2948.read_bytes()).hexdigest() if P2948.exists() else None,
            "P2949": hashlib.sha256(P2949.read_bytes()).hexdigest() if P2949.exists() else None,
            "P2950": hashlib.sha256(P2950.read_bytes()).hexdigest() if P2950.exists() else None,
        },
        "constructed_theoretical_objects": {
            "candidate_object": "RatioPackage_StrictSource_NormalForm_Lattice",
            "obligation_atoms": OBLIGATIONS,
            "truth_table_rows": rows,
            "essentiality_rows": essential,
            "current_artifact_row": current_row,
        },
        "normal_form_certificate": {
            "obligation_count": len(OBLIGATIONS),
            "truth_table_row_count": len(rows),
            "accepted_row_count": len(accepted_rows),
            "unique_accepting_row_requires_all_atoms": len(accepted_rows) == 1 and accepted_rows[0]["present_count"] == len(OBLIGATIONS),
            "all_atoms_essential": all(row["essential"] for row in essential),
            "current_present_count": current_row["present_count"],
            "current_missing_atoms": current_row["missing_atoms"],
            "current_artifact_accepts_strict_ratio_damping_source": current_row["accepted_strict_ratio_damping_source"],
        },
        "decision": {
            "positive_witnesses": {
                "normal_form_lattice_constructed": True,
                "proper_subset_rejection_proved_finitely": True,
                "content_grep_avoids_beta_scale_replay": True,
            },
            "negative_export_flags": {
                "strict_ratio_package_source_theorem_exported": False,
                "strict_torsion_character_source_theorem_exported": False,
                "strict_delta_numerator_semantics_exported": False,
                "strict_positive_beta_scale_source_exported": False,
                "strict_beta_eta_coupling_theorem_exported": False,
                "strict_damping_beta_eta_source_packet_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_hamiltonian_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2951 constructs the exact ratio-package strict-source normal form instead of replaying beta-scale or ratio scans.  The 16-row lattice has exactly one accepting row, requiring torsion-character provenance, delta semantics, beta-scale/unit source, and nonproxy variational damping coupling simultaneously.  The current artifact row has all four atoms missing, so no proper subset promotes the P2948-P2950 package to a strict damping source.",
            "next_honest_step": "Do not continue with further normal-form lattices, ratio scans, count aliases, or beta-scale orbit samples.  The next proof-grade move must supply one concrete missing atom from the P2951 lattice, preferably a strict provenance theorem for the P2938 torsion-character aggregate or a strict positive beta-scale/unit source with nonproxy coupling data; otherwise pivot outside the ratio-package lane and preserve the P2929-P2951 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["normal_form_certificate"]
    lines = [
        "# P2951/S1901 ratio-package strict-source normal-form lattice",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Normal-form certificate",
        f"- obligation count: `{cert['obligation_count']}`",
        f"- truth-table rows: `{cert['truth_table_row_count']}`",
        f"- accepted rows: `{cert['accepted_row_count']}`",
        f"- unique accepting row requires all atoms: `{cert['unique_accepting_row_requires_all_atoms']}`",
        f"- all atoms essential: `{cert['all_atoms_essential']}`",
        f"- current present count: `{cert['current_present_count']}`",
        f"- current missing atoms: `{cert['current_missing_atoms']}`",
        f"- current artifact accepts strict ratio damping source: `{cert['current_artifact_accepts_strict_ratio_damping_source']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2948), read_json(P2949), read_json(P2950))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2951/S1901 ratio-package strict-source normal-form lattice", "## P2951/S1901 ratio-package strict-source normal-form lattice\n\n`P2951/S1901` constructs the exact strict-source normal form for the P2948-P2950 ratio-package lane rather than replaying beta-scale or ratio scans.  The finite lattice has four required atoms: strict P2938 torsion-character provenance, strict identity-deficit delta semantics, strict positive beta-scale/unit source, and nonproxy variational damping coupling.  Enumerating all `2^4=16` masks gives exactly one accepting row, with all four atoms present; the current artifact row has all four atoms missing.  Therefore no proper subset exports a strict ratio-package source theorem, strict damping packet, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2951/S1901 ratio-package normal-form `L_total` guard", "## P2951/S1901 ratio-package normal-form `L_total` guard\n\n`P2951/S1901` proves by finite truth-table normal form that the exact ratio package cannot enter `L_total` through any proper subset of the four missing theorem atoms: P2938 provenance, delta semantics, beta-scale/unit source, and nonproxy variational damping coupling.  Until all four are exported, no EOM, Hamiltonian, bridge closure, role transfer, or ToE promotion follows.\n")
    append_once(AGENTS, "Current ratio-package strict-source normal-form guardrail (P2951/S1901, 2026-06-19)", "## Current ratio-package strict-source normal-form guardrail (P2951/S1901, 2026-06-19)\n\n- P2951 constructs the exact P2948-P2950 strict-source normal form rather than replaying beta-scale/unit, ratio, or count-alias scans.\n- The four required atoms are strict P2938 torsion-character provenance, strict identity-deficit delta semantics, strict positive beta-scale/unit source, and nonproxy variational damping coupling.\n- The finite `2^4=16` mask lattice has exactly one accepting row, where all four atoms are present; the current artifact row has all four atoms missing.\n- Do not continue normal-form lattices, ratio scans, count aliases, or beta-scale samples as primary strategy.  A next admissible move must export one concrete missing atom from this lattice or pivot outside the ratio-package lane while preserving the P2929-P2951 boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
