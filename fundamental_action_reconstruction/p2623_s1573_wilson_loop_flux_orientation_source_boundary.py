#!/usr/bin/env python3
from __future__ import annotations

import cmath
import hashlib
import itertools
import json
import math
import subprocess
from typing import Any

import numpy as np

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2623_s1573_wilson_loop_flux_orientation_source_boundary.json"
MD = GEN / "p2623_s1573_wilson_loop_flux_orientation_source_boundary.md"

SOURCE_FILES = {
    "P2619_SELECTOR_LATTICE": GEN / "p2619_s1569_p2618_selector_source_obligation_lattice.json",
    "P2620_TWO_OBSTRUCTION_CUT": GEN / "p2620_s1570_p2618_p2619_bridge_two_obstruction_cut.json",
    "P2621_CONDITIONAL_CHIRAL_SCHEMA": GEN / "p2621_s1571_qw636_qw1026_chiral_hopfion_selector_source_audit.json",
    "P2622_PRIOR_ART_NONPROMOTION": GEN / "p2622_s1572_qw636_qw1026_physical_rigor_nonpromotion_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "unconditional_orientation_odd_selector_source_exported",
    "unoriented_wilson_loop_promoted_to_selector_source",
    "p2620_bridge_source_cut_repaired",
    "nonlinear_damping_completion_source_exported",
    "full_bridge_revalidated",
    "role_transfer_revalidated",
    "role_bearing_ltotal_reenabled",
    "qw2191_discharged_by_this_packet",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:180]}


def semantic_rg_audit() -> dict[str, Any]:
    # Deliberately content-first patterns: no packet/ticket numbers are needed to find the research lane.
    patterns = {
        "closed_loop_flux_holonomy_content": "Wilson loop|holonomy|closed cycle|cycle product|plaquette|phase around|loop phase|flux through|Aharonov|Berry phase",
        "gauge_link_transport_content": "link phase|link variable|Peierls phase|parallel transport|gauge transform|gauge invariant|pure gauge|U_ij|connection",
        "orientation_reversal_content": "orientation reversal|cycle reversal|reverse orientation|directed loop|loop orientation|oriented cycle|unoriented cycle|sign flips",
        "selector_source_content": "orientation source|directed physical orientation|orientation-odd source|sign-sensitive physical|strict source theorem|symmetry-breaking boundary",
        "guardrail_content": "role-bearing L_total|nonlinear damping completion|bridge-source cut|QW-2191|ToE closure|strict selector closure",
    }
    return {"tool": "rg", "mode": "content-first semantic patterns, not ticket/name lookup", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def wilson_product(phases: list[float], orientation: int = 1) -> complex:
    total = sum(phases)
    return cmath.exp(1j * orientation * total)


def gauge_transform_cycle(phases: list[float], lambdas: list[float]) -> list[float]:
    n = len(phases)
    return [phases[i] + lambdas[i] - lambdas[(i + 1) % n] for i in range(n)]


def wilson_flux_certificate() -> dict[str, Any]:
    phases = [math.pi / 9.0, -math.pi / 11.0, math.pi / 7.0, math.pi / 13.0, -math.pi / 17.0]
    lambdas = [0.3, -1.2, 0.7, 2.4, -0.9]
    w_plus = wilson_product(phases, orientation=1)
    w_minus = wilson_product(phases, orientation=-1)
    transformed = gauge_transform_cycle(phases, lambdas)
    w_transformed = wilson_product(transformed, orientation=1)
    total_flux = math.remainder(sum(phases), 2.0 * math.pi)
    sign_plus = 1 if w_plus.imag > 0.0 else -1 if w_plus.imag < 0.0 else 0
    sign_minus = 1 if w_minus.imag > 0.0 else -1 if w_minus.imag < 0.0 else 0
    unoriented_class = sorted([complex(round(w_plus.real, 15), round(w_plus.imag, 15)), complex(round(w_minus.real, 15), round(w_minus.imag, 15))], key=lambda z: (z.real, z.imag))
    return {
        "analytic_identities": [
            "For U_e=exp(i a_e), W(C)=prod_e U_e=exp(i sum_e a_e).",
            "Gauge transform a_e -> a_e + lambda_tail(e)-lambda_head(e) telescopes on a closed cycle, so W(C) is gauge-invariant.",
            "Cycle reversal gives W(C^{-1})=W(C)^{-1}=conj(W(C)); therefore sign(Im W) is orientation-odd.",
            "The unoriented gauge-safe datum {W,conj(W)} or Re(W) cannot choose sign(Im W).",
        ],
        "cycle_edge_count": len(phases),
        "raw_link_phases": phases,
        "gauge_lambdas": lambdas,
        "total_flux_mod_2pi": total_flux,
        "W_oriented_real": float(w_plus.real),
        "W_oriented_imag": float(w_plus.imag),
        "W_reversed_real": float(w_minus.real),
        "W_reversed_imag": float(w_minus.imag),
        "W_after_gauge_real": float(w_transformed.real),
        "W_after_gauge_imag": float(w_transformed.imag),
        "gauge_invariance_defect_abs": abs(w_plus - w_transformed),
        "reversal_conjugacy_defect_abs": abs(w_minus - w_plus.conjugate()),
        "orientation_sign_plus": sign_plus,
        "orientation_sign_reversed": sign_minus,
        "orientation_sign_flips": sign_plus == -sign_minus and sign_plus != 0,
        "unoriented_conjugacy_class": [[z.real, z.imag] for z in unoriented_class],
        "unoriented_class_contains_both_signs": sign_plus == -sign_minus and sign_plus != 0,
        "selector_from_oriented_cycle_defined": sign_plus != 0,
        "selector_from_unoriented_cycle_defined": False,
    }


def enumerate_acceptance_rows() -> list[dict[str, Any]]:
    atoms = [
        "closed_cycle_nonzero_flux_source",
        "cycle_orientation_source",
        "gauge_safe_connection_source",
    ]
    rows = []
    for values in itertools.product([False, True], repeat=len(atoms)):
        assignment = dict(zip(atoms, values))
        orientation_atom_repaired = all(values)
        rows.append({
            "assignment": assignment,
            "orientation_odd_selector_source_repaired": orientation_atom_repaired,
            "missing": [atom for atom, value in assignment.items() if not value],
        })
    return rows


def p2620_update(orientation_repaired: bool) -> dict[str, Any]:
    rows = []
    for damping in (False, True):
        rows.append({
            "assignment": {
                "nonlinear_damping_completion_source": damping,
                "orientation_odd_selector_source_from_oriented_wilson_flux": orientation_repaired,
            },
            "bridge_source_cut_repaired": damping and orientation_repaired,
            "missing": [k for k, v in {
                "nonlinear_damping_completion_source": damping,
                "orientation_odd_selector_source_from_oriented_wilson_flux": orientation_repaired,
            }.items() if not v],
        })
    return {
        "orientation_source_exported_now": orientation_repaired,
        "rows": rows,
        "accepting_row_count": sum(1 for row in rows if row["bridge_source_cut_repaired"]),
        "verdict": "Even an oriented Wilson-flux theorem would repair only the selector atom; the current packet exports no such source and does not repair nonlinear damping.",
    }


def build_payload() -> dict[str, Any]:
    artifacts = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    flux = wilson_flux_certificate()
    acceptance_rows = enumerate_acceptance_rows()
    current_orientation_repaired = False
    theorem_export = {
        "theorem_name": "P2623 oriented Wilson-flux selector boundary theorem",
        "claim": (
            "A gauge-invariant Wilson loop can supply an orientation-odd sign only after three typed inputs are present: a gauge-safe connection on a closed cycle, nonzero flux with Im W != 0, and an independently sourced orientation of that cycle.  The gauge-invariant but unoriented datum {W,W^{-1}} is sign-blind, so Wilson/holonomy content alone does not repair the P2620 orientation atom."
        ),
        "positive_content": [
            "Closed-cycle Wilson products are gauge-invariant by telescoping link gauge shifts.",
            "Reversing the oriented cycle conjugates the Wilson loop and flips sign(Im W) when the flux is not 0 or pi.",
            "If a future strict theorem supplies gauge-safe connection, nonzero flux, and cycle orientation, sigma=sign(Im W) is a mathematically valid orientation-odd selector candidate.",
        ],
        "obstructions": [
            "Without a cycle-orientation source the physical datum is the conjugacy pair {W,W^{-1}}, which contains both signs.",
            "A Wilson loop can be a gauge-invariant flux diagnostic, but it cannot create the missing orientation by itself.",
            "No nonlinear damping completion source is supplied by this selector-side analysis.",
        ],
        "next_admissible_research_targets": [
            "search for an internal strict source of oriented closed cycles rather than another holonomy value",
            "derive a gauge-safe connection/parallel transport from strict nadsoliton field content",
            "keep nonlinear damping completion as a separate bridge atom",
        ],
        "not_licensed": [
            "unconditional orientation_odd_selector_source",
            "promotion of unoriented Wilson-loop data to strict selector source",
            "P2620 bridge-source cut repair",
            "role-bearing L_total",
            "QW-2191 discharge",
            "ToE closure",
        ],
    }
    payload = {
        "packet_id": "P2623",
        "slice_id": "S1573",
        "status": "P2623_WILSON_FLUX_ORIENTATION_SOURCE_BOUNDARY_NO_SELECTOR_SOURCE_NO_BRIDGE_NO_LTOTAL_NO_QW2191_NO_TOE",
        "source_artifacts": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "inherited_artifact_status": {name: artifact.get("status", artifact.get("packet_id", "UNKNOWN")) for name, artifact in artifacts.items()},
        "wilson_flux_certificate": flux,
        "orientation_source_acceptance_lattice": {
            "atoms": ["closed_cycle_nonzero_flux_source", "cycle_orientation_source", "gauge_safe_connection_source"],
            "rows": acceptance_rows,
            "accepting_row_count": sum(1 for row in acceptance_rows if row["orientation_odd_selector_source_repaired"]),
            "current_orientation_repaired": current_orientation_repaired,
        },
        "p2620_update": p2620_update(current_orientation_repaired),
        "theorem_export": theorem_export,
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    payload["certificate_hash"] = sha256_json({k: v for k, v in payload.items() if k != "certificate_hash"})
    return payload


def write_markdown(payload: dict[str, Any]) -> None:
    theorem = payload["theorem_export"]
    flux = payload["wilson_flux_certificate"]
    lines = [
        "# P2623/S1573 Wilson-loop flux orientation-source boundary",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication grep audit",
        "",
        f"Mode: `{payload['semantic_rg_antiduplication_audit']['mode']}`.",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits; samples retained in JSON certificate.")
    lines.extend([
        "",
        "## Theorem export",
        "",
        f"**Claim.** {theorem['claim']}",
        "",
        "Positive retained content:",
    ])
    for item in theorem["positive_content"]:
        lines.append(f"- {item}")
    lines.append("")
    lines.append("Obstructions:")
    for item in theorem["obstructions"]:
        lines.append(f"- {item}")
    lines.extend([
        "",
        "## Computational certificate",
        "",
        f"- Gauge-invariance defect: `{flux['gauge_invariance_defect_abs']}`.",
        f"- Reversal conjugacy defect: `{flux['reversal_conjugacy_defect_abs']}`.",
        f"- Oriented sign flips under reversal: `{flux['orientation_sign_flips']}`.",
        f"- Unoriented conjugacy class contains both signs: `{flux['unoriented_class_contains_both_signs']}`.",
        f"- Orientation-source lattice accepting rows: `{payload['orientation_source_acceptance_lattice']['accepting_row_count']}` of 8.",
        f"- P2620 accepting rows now: `{payload['p2620_update']['accepting_row_count']}`.",
        "",
        "## Next admissible targets",
    ])
    for item in theorem["next_admissible_research_targets"]:
        lines.append(f"- {item}")
    lines.extend(["", "Not licensed:"])
    for item in theorem["not_licensed"]:
        lines.append(f"- {item}")
    lines.extend(["", f"Certificate hash: `{payload['certificate_hash']}`", ""])
    MD.write_text("\n".join(lines), encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2623/S1573 Wilson-loop flux orientation-source boundary

`P2623/S1573` performs the next content-first selector audit: Wilson/holonomy data are useful only after the missing typed sources are separated.  A closed-cycle Wilson product is gauge-invariant and an oriented nonzero flux has `sign(Im W)` that flips under cycle reversal, but the unoriented gauge-safe datum is the conjugacy pair `{W,W^{-1}}` and is sign-blind.  Therefore Wilson-loop/flux content alone does not export `orientation_odd_selector_source`; it would require a gauge-safe connection source, nonzero flux source, and independent cycle-orientation source, and it still would not supply nonlinear damping completion, role-bearing `L_total`, QW-2191 discharge, or ToE closure.
"""
    lag_section = """
## P2623/S1573 Wilson-loop orientation-source Ltotal guard

`P2623/S1573` keeps role-bearing `L_total` closed.  The Wilson-loop calculation identifies a possible future selector shape, `sigma=sign(Im W)`, but only under a real oriented-cycle/connection/nonzero-flux source theorem.  Unoriented holonomy data remain sign-blind, and the independent nonlinear damping completion atom remains open.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2623/S1573 Wilson-loop flux orientation-source boundary", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2623/S1573 Wilson-loop orientation-source Ltotal guard", lag_section)


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps({
        "packet_id": payload["packet_id"],
        "status": payload["status"],
        "gauge_invariance_defect": payload["wilson_flux_certificate"]["gauge_invariance_defect_abs"],
        "reversal_conjugacy_defect": payload["wilson_flux_certificate"]["reversal_conjugacy_defect_abs"],
        "orientation_sign_flips": payload["wilson_flux_certificate"]["orientation_sign_flips"],
        "orientation_lattice_accepting_rows": payload["orientation_source_acceptance_lattice"]["accepting_row_count"],
        "p2620_accepting_rows_now": payload["p2620_update"]["accepting_row_count"],
        "certificate_hash": payload["certificate_hash"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
