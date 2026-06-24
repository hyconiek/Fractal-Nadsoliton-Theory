#!/usr/bin/env python3
"""P3070/S2020: canonical length/time unit-provider obstruction.

P3069 left one honest blocker: a canonical length/time unit provider that
converts a nonconventional nadsoliton scalar into a unit-bearing coordinate.
This step attacks exactly that unit-source atom.  It does not replay selector
closure and does not promote chart rank, spacetime embedding, or observed light.

The constructed object is a finite provider matrix.  Candidate nadsoliton
scalars are separated from candidate unit maps and from the strict acceptance
requirements: intrinsic scalar, canonical positive scale, nonconventional unit
law, length/time dual availability, and coordinate coupling back to P3069.
"""
from __future__ import annotations

import hashlib, json, subprocess
from fractions import Fraction
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3069_s2019_coordinate_pair_source_rank_obstruction import OUT as P3069

OUT = GEN / "p3070_s2020_canonical_length_time_unit_provider_obstruction.json"
MD = GEN / "p3070_s2020_canonical_length_time_unit_provider_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "canonical_length_time_unit_provider": r"canonical.*(length|time).*unit|length/time unit|unit[- ]source|unit provider",
    "nadsoliton_scalar_to_unit": r"nadsoliton.*scalar|nonconventional.*scalar|scalar.*unit|scale[- ]control",
    "coordinate_unit_coupling": r"unit-bearing.*coordinate|t_sigma|x_sigma|coordinate.*unit|speed[- ]of[- ]light",
    "observed_light_boundaries": r"observed light|gauge photon|Lorentz closure|empirical physics|unit-bearing action",
}

SCALARS = (
    {"id": "sigma_branch_sign", "value": Fraction(1), "intrinsic": True, "positive": False, "nonconventional": False, "blocked_by": "sign branch is axiom-conditioned and not a positive unit scale"},
    {"id": "z12_orbit_cardinality", "value": Fraction(12), "intrinsic": True, "positive": True, "nonconventional": False, "blocked_by": "cardinality is dimensionless counting data, not a length/time unit law"},
    {"id": "alpha_geo_shape", "value": "4*ln(2)", "intrinsic": True, "positive": True, "nonconventional": False, "blocked_by": "shape normalization lacks role-safe unit absorption/source theorem"},
    {"id": "strict_gate_beta", "value": Fraction(1), "intrinsic": False, "positive": True, "nonconventional": False, "blocked_by": "later gate control is not a canonical ontology-side length/time unit"},
    {"id": "entropy_one_bit", "value": Fraction(1), "intrinsic": False, "positive": True, "nonconventional": False, "blocked_by": "missing intrinsic reference cell and bit-to-length/action map"},
    {"id": "formal_unit_slot", "value": None, "intrinsic": False, "positive": False, "nonconventional": False, "blocked_by": "placeholder has no exported scalar value"},
)

UNIT_MAPS = (
    {"id": "identity_dimensionless_scale", "length_unit": False, "time_unit": False, "nonconventional_unit_law": False, "coordinate_coupling": False},
    {"id": "declare_length_l0", "length_unit": True, "time_unit": False, "nonconventional_unit_law": False, "coordinate_coupling": False},
    {"id": "declare_time_t0", "length_unit": False, "time_unit": True, "nonconventional_unit_law": False, "coordinate_coupling": False},
    {"id": "declare_c_equals_one_bridge", "length_unit": True, "time_unit": True, "nonconventional_unit_law": False, "coordinate_coupling": False},
    {"id": "p3069_coordinate_rerun_hook", "length_unit": False, "time_unit": False, "nonconventional_unit_law": False, "coordinate_coupling": True},
)


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:10]})
    return rows


def provider_matrix() -> list[dict[str, Any]]:
    rows = []
    for scalar, unit_map in product(SCALARS, UNIT_MAPS):
        accepted = all([
            scalar["intrinsic"], scalar["positive"], scalar["nonconventional"],
            unit_map["length_unit"], unit_map["time_unit"], unit_map["nonconventional_unit_law"], unit_map["coordinate_coupling"],
        ])
        rows.append({
            "scalar_id": scalar["id"],
            "unit_map_id": unit_map["id"],
            "intrinsic_scalar": scalar["intrinsic"],
            "positive_scale": scalar["positive"],
            "nonconventional_scalar_source": scalar["nonconventional"],
            "length_unit_available": unit_map["length_unit"],
            "time_unit_available": unit_map["time_unit"],
            "nonconventional_unit_law": unit_map["nonconventional_unit_law"],
            "coordinate_coupling_to_p3069": unit_map["coordinate_coupling"],
            "accepted_canonical_length_time_unit_provider": accepted,
            "blocked_by": "; ".join(filter(None, [
                None if scalar["intrinsic"] else "scalar not intrinsic to current nadsoliton source layer",
                None if scalar["positive"] else "no positive scale value",
                None if scalar["nonconventional"] else scalar["blocked_by"],
                None if unit_map["length_unit"] else "no length unit",
                None if unit_map["time_unit"] else "no time unit",
                None if unit_map["nonconventional_unit_law"] else "unit map is declaration/convention, not strict source law",
                None if unit_map["coordinate_coupling"] else "no coupling to P3069 t_sigma/x_sigma coordinate source",
            ])),
        })
    return rows


def build_payload() -> dict[str, Any]:
    p3069 = read_json(P3069)
    grep_rows = content_grep()
    rows = provider_matrix()
    accepted = [r for r in rows if r["accepted_canonical_length_time_unit_provider"]]
    proof_obligations = [
        {"obligation": "content_first_grep_before_unit_provider_audit", "satisfied": True, "detail": "searched by unit-provider, scalar-to-unit, coordinate-unit coupling, and observed-light boundary content"},
        {"obligation": "construct_candidate_scalar_and_unit_map_space", "satisfied": True, "detail": "six scalar candidates crossed with five candidate unit maps"},
        {"obligation": "finite_unit_provider_matrix", "satisfied": True, "detail": "30 rows audited against intrinsic-positive-nonconventional-unit-coordinate criteria"},
        {"obligation": "separate_dimensionless_scale_from_unit_bearing_coordinate", "satisfied": True, "detail": "positive dimensionless scalars are not accepted without nonconventional length/time unit law and coordinate coupling"},
        {"obligation": "export_canonical_length_time_unit_provider", "satisfied": False, "detail": "no current scalar-map pair satisfies all strict unit-source criteria"},
    ]
    return {
        "status": "P3070_CANONICAL_LENGTH_TIME_UNIT_PROVIDER_OBSTRUCTION_NO_EXPORT",
        "input_hashes": {"P3069": hashlib.sha256(P3069.read_bytes()).hexdigest() if P3069.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "canonical_unit_provider_template": {
                "object": "CanonicalLengthTimeUnitProviderTemplate",
                "formal_shape": "U(s) -> (ell_s, tau_s, c_s=ell_s/tau_s) with coupling to P3069 t_sigma/x_sigma",
                "acceptance_criteria": ["intrinsic_scalar", "positive_scale", "nonconventional_scalar_source", "length_unit", "time_unit", "nonconventional_unit_law", "coordinate_coupling_to_p3069"],
            },
            "candidate_scalars": [{**s, "value": str(s["value"]) if s["value"] is not None else None} for s in SCALARS],
            "candidate_unit_maps": list(UNIT_MAPS),
            "unit_provider_obstruction_matrix": rows,
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(r["hit_count"] for r in grep_rows),
            "candidate_scalars": len(SCALARS),
            "candidate_unit_maps": len(UNIT_MAPS),
            "provider_matrix_rows": len(rows),
            "positive_scalar_candidates": sum(1 for s in SCALARS if s["positive"]),
            "intrinsic_scalar_candidates": sum(1 for s in SCALARS if s["intrinsic"]),
            "nonconventional_scalar_candidates": sum(1 for s in SCALARS if s["nonconventional"]),
            "unit_maps_with_length_and_time": sum(1 for u in UNIT_MAPS if u["length_unit"] and u["time_unit"]),
            "unit_maps_with_coordinate_coupling": sum(1 for u in UNIT_MAPS if u["coordinate_coupling"]),
            "nonconventional_unit_laws": sum(1 for u in UNIT_MAPS if u["nonconventional_unit_law"]),
            "accepted_unit_provider_rows": len(accepted),
            "p3069_accepted_coordinate_pair_source_rows": p3069.get("finite_certificate", {}).get("accepted_coordinate_pair_source_rows"),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(1 for o in proof_obligations if o["satisfied"]),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3070 attacks the exact P3069 unit-source atom by constructing a CanonicalLengthTimeUnitProviderTemplate and a 6 x 5 finite scalar/unit-map matrix.  Current artifacts contain positive dimensionless scalars and conventional declarations, but 0 nonconventional scalar sources, 0 nonconventional unit laws, and 0 accepted provider rows that simultaneously export length, time, and coordinate coupling.  Thus P3069 cannot yet be rerun as a unit-bearing coordinate-source theorem.",
            "negative_export_flags": {k: False for k in ["canonical_length_time_unit_provider_exported", "unit_bearing_coordinate_exported", "coordinate_pair_source_theorem_exported", "strict_spacetime_embedding_exported", "observed_light_exported", "gauge_photon_sector_exported", "unit_bearing_action_exported", "variational_EOM_exported", "empirical_physics_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"unit_provider_template_constructed": True, "finite_scalar_unit_matrix_executed": True, "dimensionless_vs_unit_bearing_boundary_separated": True},
            "next_honest_step": "Do not declare units by convention or promote dimensionless scalars to spacetime.  Since the P3069 unit-source atom is obstructed on current artifacts, pivot to the P3066 sigma-invariant scalar conservation/scale-control row: construct one finite conservation or bounded-scale-control theorem for a sigma-even nadsoliton scalar, with no observed-light, selector, L_total, bridge/role-transfer, or ToE promotion.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3070/S2020 canonical length/time unit-provider obstruction", "", f"Status: `{payload['status']}`", "",
        "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- candidate scalars: `{c['candidate_scalars']}`",
        f"- candidate unit maps: `{c['candidate_unit_maps']}`",
        f"- provider matrix rows: `{c['provider_matrix_rows']}`",
        f"- positive scalar candidates: `{c['positive_scalar_candidates']}`",
        f"- intrinsic scalar candidates: `{c['intrinsic_scalar_candidates']}`",
        f"- nonconventional scalar candidates: `{c['nonconventional_scalar_candidates']}`",
        f"- unit maps with length and time: `{c['unit_maps_with_length_and_time']}`",
        f"- unit maps with coordinate coupling: `{c['unit_maps_with_coordinate_coupling']}`",
        f"- nonconventional unit laws: `{c['nonconventional_unit_laws']}`",
        f"- accepted unit provider rows: `{c['accepted_unit_provider_rows']}`",
        f"- P3069 accepted coordinate-pair source rows: `{c['p3069_accepted_coordinate_pair_source_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3070/S2020 canonical length/time unit-provider obstruction", "## P3070/S2020 canonical length/time unit-provider obstruction\n\n`P3070/S2020` attacks one P3069 atom: a canonical length/time unit provider that would convert a nonconventional nadsoliton scalar into a unit-bearing coordinate.  It constructs a `CanonicalLengthTimeUnitProviderTemplate` and audits `30` scalar/unit-map rows.  Current artifacts include positive dimensionless candidates and conventional declarations, but `0` nonconventional scalar candidates, `0` nonconventional unit laws, and `0` accepted unit-provider rows.  Therefore no unit-bearing coordinate, strict spacetime embedding, observed-light closure, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3070/S2020 unit declarations are not action units", "## P3070/S2020 unit declarations are not action units\n\n`P3070/S2020` blocks the move from dimensionless nadsoliton scalars or declared `c=1` conventions to a unit-bearing coordinate/action source.  The audit exports no length/time unit law, no light-sector action density, no EOM/Hamiltonian, and no empirical constant map.\n")
    append_once(AGENTS, "Current canonical length/time unit-provider obstruction guardrail (P3070/S2020, 2026-06-24)", "## Current canonical length/time unit-provider obstruction guardrail (P3070/S2020, 2026-06-24)\n\n- P3070 attacks exactly one P3069 atom: a canonical length/time unit provider for converting a nonconventional nadsoliton scalar into a unit-bearing coordinate.\n- The finite scalar/unit-map audit has `6` scalar candidates, `5` unit-map candidates, and `30` provider rows; positive dimensionless candidates exist, but there are `0` nonconventional scalar candidates, `0` nonconventional unit laws, and `0` accepted canonical length/time unit-provider rows.\n- Do not promote dimensionless scalars or conventional unit declarations to unit-bearing coordinates, strict spacetime embedding, observed light, gauge photons, unit-bearing action, variational EOM, empirical physics, `QW-2191` discharge, strict selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is to pivot to the P3066 sigma-invariant scalar conservation/scale-control row: construct one finite conservation or bounded-scale-control theorem for a sigma-even nadsoliton scalar, without observed-light or closure promotion.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
