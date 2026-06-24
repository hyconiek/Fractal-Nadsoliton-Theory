#!/usr/bin/env python3
"""P3072/S2022: discrete continuity/Noether-current interface matrix.

P3071 exported scoped dimensionless sigma-even scalar conservation rows and
recommended exactly one transition-interface theorem toward dynamics.  P3072
constructs that interface on the finite Z12 cycle: scalar densities, edge
currents, an exact incidence/divergence operator, and a continuity equation

    delta_t rho + div J = 0.

The result is deliberately bounded.  Static scalar conservation admits the zero
current without extra premises.  Nonzero cycle currents solve the divergence
equation, but their sign/orientation is premise-based on the current artifacts;
gradient/shell currents fail divergence.  Therefore P3072 exports an exact
interface object and a no-go for nontrivial non-premise Noether-current dynamics,
not spacetime, units, observed light, gauge photons, EOM, L_total, bridge/role
transfer, selector closure, or ToE closure.
"""
from __future__ import annotations

import hashlib, json, subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3071_s2021_sigma_invariant_scalar_conservation_scale_control import OUT as P3071, profile_values, scalar_summary, transform_values

OUT = GEN / "p3072_s2022_discrete_continuity_noether_current_matrix.json"
MD = GEN / "p3072_s2022_discrete_continuity_noether_current_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "discrete_continuity_noether_interface": r"discrete continuity|Noether[- ]current|continuity equation|current matrix",
    "sigma_even_scalar_input": r"sigma-even.*scalar|sigma[- ]invariant.*scalar|bounded.*scale-control",
    "z12_incidence_divergence": r"Z12.*incidence|divergence|cycle current|edge current",
    "no_dynamics_promotion_boundary": r"observed light|gauge photon|unit-bearing action|variational EOM|empirical physics|L_total|ToE|selector closure",
}

SIGMAS = (-1, 1)
Z12 = tuple(range(12))
TRANSFORMS = tuple((kind, k) for kind in ("rotation", "reflection") for k in Z12)
ACCEPTED_P3071_PROFILES = (
    "constant_cardinality_density",
    "even_distance_quadratic_density",
    "even_distance_shell_indicator_density",
)
CURRENT_TEMPLATES = (
    {"id": "zero_current", "orientation_premise": False, "nontrivial": False},
    {"id": "oriented_constant_cycle_current", "orientation_premise": True, "nontrivial": True},
    {"id": "sigma_oriented_constant_cycle_current", "orientation_premise": True, "nontrivial": True},
    {"id": "forward_gradient_current", "orientation_premise": False, "nontrivial": True},
    {"id": "alternating_shell_current", "orientation_premise": False, "nontrivial": True},
)


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:10]})
    return rows


def current_values(template_id: str, rho: tuple[Fraction, ...], sigma: int) -> tuple[Fraction, ...]:
    if template_id == "zero_current":
        return tuple(Fraction(0) for _ in Z12)
    if template_id == "oriented_constant_cycle_current":
        return tuple(Fraction(1) for _ in Z12)
    if template_id == "sigma_oriented_constant_cycle_current":
        return tuple(Fraction(sigma) for _ in Z12)
    if template_id == "forward_gradient_current":
        return tuple(rho[(n + 1) % 12] - rho[n] for n in Z12)
    if template_id == "alternating_shell_current":
        return tuple(Fraction(1 if n % 2 == 0 else -1) for n in Z12)
    raise ValueError(template_id)


def divergence(edge_current: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    # Edge n is n -> n+1.  Net outflow at vertex n is J_n - J_{n-1}.
    return tuple(edge_current[n] - edge_current[(n - 1) % 12] for n in Z12)


def support_size(values: tuple[Fraction, ...]) -> int:
    return sum(1 for value in values if value != 0)


def matrix_rows() -> list[dict[str, Any]]:
    rows = []
    for profile_id in ACCEPTED_P3071_PROFILES:
        for sigma in SIGMAS:
            base_rho = profile_values(profile_id, sigma)
            for kind, k in TRANSFORMS:
                rho = transform_values(base_rho, kind, k)
                rho_summary = scalar_summary(rho)
                for template in CURRENT_TEMPLATES:
                    current = current_values(template["id"], rho, sigma)
                    div = divergence(current)
                    continuity_residual = div  # static delta_t rho = 0 at this interface stage.
                    divergence_zero = all(v == 0 for v in div)
                    accepted_interface_row = divergence_zero and not template["orientation_premise"] and not template["nontrivial"]
                    accepted_nontrivial_noether_row = divergence_zero and not template["orientation_premise"] and template["nontrivial"] and support_size(current) > 0
                    rows.append({
                        "profile_id": profile_id,
                        "sigma": sigma,
                        "transform": {"kind": kind, "k": k},
                        "current_template": template["id"],
                        "rho_summary": rho_summary,
                        "current_support": support_size(current),
                        "divergence_support": support_size(div),
                        "continuity_residual_support": support_size(continuity_residual),
                        "divergence_zero": divergence_zero,
                        "orientation_premise_required": template["orientation_premise"],
                        "nontrivial_current": template["nontrivial"],
                        "accepted_premise_free_static_continuity_row": accepted_interface_row,
                        "accepted_nontrivial_noether_current_row": accepted_nontrivial_noether_row,
                        "blocked_by": "" if accepted_interface_row or accepted_nontrivial_noether_row else "; ".join(filter(None, [
                            None if divergence_zero else "nonzero divergence/continuity residual",
                            "orientation premise required" if template["orientation_premise"] else None,
                            "trivial zero-current only" if not template["nontrivial"] else None,
                            "no non-premise dynamical update law" if divergence_zero and template["nontrivial"] else None,
                        ])),
                    })
    return rows


def build_payload() -> dict[str, Any]:
    p3071 = read_json(P3071)
    rows = matrix_rows()
    grep_rows = content_grep()
    accepted_static = [r for r in rows if r["accepted_premise_free_static_continuity_row"]]
    accepted_nontrivial = [r for r in rows if r["accepted_nontrivial_noether_current_row"]]
    divergence_zero = [r for r in rows if r["divergence_zero"]]
    proof_obligations = [
        {"obligation": "content_first_grep_before_current_interface", "satisfied": True, "detail": "searched by continuity/Noether, sigma scalar, Z12 divergence, and no-promotion content"},
        {"obligation": "construct_z12_incidence_divergence_operator", "satisfied": True, "detail": "edge current J_n on C12 with div(J)_n = J_n - J_{n-1}"},
        {"obligation": "test_p3071_profiles_against_current_templates", "satisfied": True, "detail": "3 profiles x 2 sigma branches x 24 D12 transforms x 5 current templates = 720 exact rows"},
        {"obligation": "export_premise_free_static_continuity_interface", "satisfied": bool(accepted_static), "detail": f"{len(accepted_static)} zero-current rows satisfy static continuity without orientation premises"},
        {"obligation": "export_nontrivial_noether_current_dynamics", "satisfied": False, "detail": "nonzero divergence-free cycle currents require an orientation/sign premise; non-oriented gradient/shell currents have nonzero divergence"},
    ]
    return {
        "status": "P3072_DISCRETE_CONTINUITY_NOETHER_CURRENT_INTERFACE_BOUNDED_NO_GO",
        "input_hashes": {"P3071": hashlib.sha256(P3071.read_bytes()).hexdigest() if P3071.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "discrete_continuity_noether_current_template": {
                "object": "Z12DiscreteContinuityNoetherCurrentInterfaceTemplate",
                "domain": "P3071 accepted sigma-even scalar profiles over T_sigma x Z12",
                "vertices": 12,
                "directed_edges": 12,
                "incidence_rule": "div(J)_n = J_n - J_{n-1}",
                "continuity_equation": "delta_t rho + div J = 0 with delta_t rho fixed to 0 at this interface stage",
            },
            "current_templates": list(CURRENT_TEMPLATES),
            "continuity_current_matrix": rows,
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(r["hit_count"] for r in grep_rows),
            "p3071_accepted_profiles": p3071.get("finite_certificate", {}).get("accepted_profile_count"),
            "profiles_tested": len(ACCEPTED_P3071_PROFILES),
            "sigma_branches": len(SIGMAS),
            "d12_transforms": len(TRANSFORMS),
            "current_templates": len(CURRENT_TEMPLATES),
            "continuity_matrix_rows": len(rows),
            "divergence_zero_rows": len(divergence_zero),
            "orientation_premise_rows": sum(1 for r in rows if r["orientation_premise_required"]),
            "accepted_premise_free_static_continuity_rows": len(accepted_static),
            "accepted_nontrivial_noether_current_rows": len(accepted_nontrivial),
            "nonzero_divergence_rows": sum(1 for r in rows if not r["divergence_zero"]),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(1 for o in proof_obligations if o["satisfied"]),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3072 constructs the requested transition-interface object from P3071 conserved scalars to a candidate dynamics layer: an exact Z12 incidence/divergence operator and a finite continuity-current matrix.  The only premise-free accepted rows are static zero-current rows.  Nonzero cycle currents are divergence-free but orientation/sign-premise based, while gradient and alternating-shell current templates have nonzero divergence.  Thus no nontrivial non-premise Noether-current dynamics is exported on current artifacts.",
            "negative_export_flags": {k: False for k in ["nontrivial_noether_current_exported", "dynamical_update_law_exported", "canonical_length_time_unit_provider_exported", "unit_bearing_coordinate_exported", "strict_spacetime_embedding_exported", "observed_light_exported", "gauge_photon_sector_exported", "unit_bearing_action_exported", "variational_EOM_exported", "empirical_physics_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"z12_incidence_divergence_operator_exported": True, "premise_free_static_continuity_interface_exported": True, "finite_current_matrix_executed": True},
            "next_honest_step": "Do not replay selector or promote the zero-current interface to physics.  The next proof-grade move is one bounded renormalization/scale-flow obstruction table for the P3071 sigma-even scalar summaries: test whether any intrinsic, premise-free scale-flow operator preserves the accepted summaries while producing a nonzero bounded flow.  If that also fails, preserve the bounded no-dynamics certificate until a new strict action/EOM/unit provider is introduced.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3072/S2022 discrete continuity/Noether-current interface matrix", "", f"Status: `{payload['status']}`", "",
        "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- P3071 accepted profiles: `{c['p3071_accepted_profiles']}`",
        f"- profiles tested: `{c['profiles_tested']}`",
        f"- sigma branches: `{c['sigma_branches']}`",
        f"- D12 transforms: `{c['d12_transforms']}`",
        f"- current templates: `{c['current_templates']}`",
        f"- continuity matrix rows: `{c['continuity_matrix_rows']}`",
        f"- divergence-zero rows: `{c['divergence_zero_rows']}`",
        f"- orientation-premise rows: `{c['orientation_premise_rows']}`",
        f"- accepted premise-free static continuity rows: `{c['accepted_premise_free_static_continuity_rows']}`",
        f"- accepted nontrivial Noether-current rows: `{c['accepted_nontrivial_noether_current_rows']}`",
        f"- nonzero-divergence rows: `{c['nonzero_divergence_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3072/S2022 discrete continuity/Noether-current interface", "## P3072/S2022 discrete continuity/Noether-current interface\n\n`P3072/S2022` constructs `Z12DiscreteContinuityNoetherCurrentInterfaceTemplate`: an exact `Z12` incidence/divergence operator `div(J)_n = J_n - J_{n-1}` and a finite continuity matrix over the three `P3071` accepted sigma-even scalar profiles.  It audits `3` profiles across `2` sigma branches, `24` `D12` transforms, and `5` current templates, for `720` exact rows.  Only zero-current static continuity is premise-free; nonzero cycle currents require an orientation/sign premise and gradient/shell currents have nonzero divergence.  The result exports no nontrivial Noether current, unit-bearing action, EOM, observed-light/gauge sector, `L_total`, bridge/role-transfer, selector closure, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3072/S2022 continuity interface is not an EOM", "## P3072/S2022 continuity interface is not an EOM\n\n`P3072/S2022` adds an exact finite continuity-current interface for the `P3071` dimensionless scalar profiles.  The interface proves that static zero current is premise-free, while nonzero divergence-free cycle currents still require an orientation/sign premise.  This is not a Noether theorem, action density, Hamiltonian, variational EOM, gauge field, observed spacetime dynamics, or empirical physics map.\n")
    append_once(AGENTS, "Current discrete continuity/Noether-current interface guardrail (P3072/S2022, 2026-06-24)", "## Current discrete continuity/Noether-current interface guardrail (P3072/S2022, 2026-06-24)\n\n- P3072 follows the P3071 recommendation and constructs one transition-interface theorem: an exact `Z12` incidence/divergence operator and finite continuity-current matrix for the accepted sigma-even scalar profiles.\n- The matrix has `3` profiles, `2` sigma branches, `24` `D12` transforms, `5` current templates, and `720` exact rows.  Static zero current is premise-free; nonzero cycle currents are divergence-free but orientation/sign-premise based; gradient/shell current templates have nonzero divergence.\n- Do not promote P3072 to a nontrivial Noether current, dynamical update law, canonical unit provider, unit-bearing coordinates, strict spacetime embedding, observed light, gauge photons, unit-bearing action, variational EOM, empirical physics, `QW-2191` discharge, strict selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is one bounded renormalization/scale-flow obstruction table for the P3071 scalar summaries, unless a new strict action/EOM/unit provider is introduced.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
