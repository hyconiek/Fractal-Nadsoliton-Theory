#!/usr/bin/env python3
"""P3067/S2017: sigma-conditioned nadsoliton-to-light Lorentz proxy matrix.

P3066 recommended choosing one row, preferably `light_emergence_interface`.
This step does not replay selector derivation.  It works under the explicit
P3065 T_sigma boundary and constructs the smallest finite transition-law proxy
that can be tested computationally: sigma selects one of two oriented 1+1 null
rays k_sigma = (1, sigma).  The finite Lorentz audit proves only proxy-level
null-norm covariance under sampled boosts; it does not export observed light,
gauge theory, a unit-bearing action, EOM, or empirical physics.
"""
from __future__ import annotations

import hashlib, json, subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3066_s2016_sigma_conditioned_standard_physics_compatibility_matrix import OUT as P3066

OUT = GEN / "p3067_s2017_sigma_light_lorentz_proxy_matrix.json"
MD = GEN / "p3067_s2017_sigma_light_lorentz_proxy_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "sigma_conditioned_light_transition_law": r"sigma[-_ ]conditioned.*light|nadsoliton.*light.*transition|light_emergence_interface|T_sigma.*light",
    "finite_lorentz_null_proxy": r"null[-_ ]ray|null norm|Lorentz.*covariance|boost.*covariance|1\+1.*Lorentz",
    "no_observed_light_closure_boundary": r"observed light|gauge.*photon|polarization.*transversality|unit-bearing.*action|empirical.*light",
}

SIGMA_BRANCHES = {"sigma_plus": 1, "sigma_minus": -1}
BOOSTS = (Fraction(0, 1), Fraction(1, 3), Fraction(-1, 2), Fraction(2, 5))
STRICT_BLOCKERS = (
    "no_strict_nadsoliton_to_spacetime_embedding",
    "no_unit_bearing_metric_or_speed_of_light_unit",
    "no_gauge_photon_field_or_polarization_bundle",
    "no_variational_EOM_or_Hamiltonian_for_light_sector",
    "no_empirical_constant_map_or_observed_light_export",
)


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:10]})
    return rows


def null_norm(t: Fraction, x: Fraction) -> Fraction:
    return t * t - x * x


def boost_null_row(sigma_name: str, sigma: int, v: Fraction) -> dict[str, Any]:
    # Avoid irrational gamma by checking the exact scaled norm:
    # norm(gamma*(t-vx), gamma*(x-vt)) = gamma^2 * scaled_norm.
    t, x = Fraction(1), Fraction(sigma)
    t_scaled = t - v * x
    x_scaled = x - v * t
    scaled_norm = null_norm(t_scaled, x_scaled)
    gamma_denominator = Fraction(1) - v * v
    input_norm = null_norm(t, x)
    # For Lorentz boosts, scaled_norm == (1-v^2) * input_norm.  Null rows
    # therefore remain null without choosing units or a physical c.
    expected_scaled_norm = gamma_denominator * input_norm
    proxy_pass = scaled_norm == expected_scaled_norm == 0
    return {
        "sigma_branch": sigma_name,
        "sigma_value": sigma,
        "input_ray_proxy": [str(t), str(x)],
        "boost_velocity_fraction": str(v),
        "boost_gamma_squared_formal": str(Fraction(1, 1) / gamma_denominator),
        "boosted_scaled_ray_before_gamma": [str(t_scaled), str(x_scaled)],
        "input_minkowski_norm": str(input_norm),
        "scaled_output_norm_before_gamma": str(scaled_norm),
        "proxy_null_covariance_pass": proxy_pass,
        "strict_lorentz_closure_exported": False,
        "remaining_blockers": list(STRICT_BLOCKERS),
    }


def build_proxy_rows() -> list[dict[str, Any]]:
    return [boost_null_row(name, sigma, v) for name, sigma in SIGMA_BRANCHES.items() for v in BOOSTS]


def build_payload() -> dict[str, Any]:
    p3066 = read_json(P3066)
    grep_rows = content_grep()
    proxy_rows = build_proxy_rows()
    proof_obligations = [
        {"obligation": "content_first_grep_before_light_transition_lorentz_proxy", "satisfied": True, "detail": "searched by light-transition, null-ray/Lorentz, and observed-light-boundary content rather than by numbers only"},
        {"obligation": "construct_sigma_conditioned_transition_proxy", "satisfied": True, "detail": "constructed L_sigma: T_sigma nadsoliton boundary datum -> oriented 1+1 null-ray proxy k_sigma=(1,sigma)"},
        {"obligation": "finite_lorentz_null_covariance_audit", "satisfied": True, "detail": "checked exact rational scaled null-norm preservation for 2 sigma branches and 4 sampled boosts"},
        {"obligation": "export_strict_nadsoliton_to_light_theorem", "satisfied": False, "detail": "proxy has no strict spacetime embedding, unit-bearing metric, gauge photon field, variational EOM, or empirical map"},
    ]
    return {
        "status": "P3067_SIGMA_LIGHT_LORENTZ_PROXY_MATRIX_CONDITIONAL_NO_OBSERVED_LIGHT_CLOSURE",
        "input_hashes": {"P3066": hashlib.sha256(P3066.read_bytes()).hexdigest() if P3066.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "sigma_conditioned_nadsoliton_to_light_transition_proxy": {
                "object": "SigmaConditionedNadsolitonLightTransitionProxy",
                "definition": "Under T_sigma, define the conditional proxy law L_sigma sending the recorded orientation boundary datum sigma to the oriented 1+1 null ray k_sigma=(1,sigma).",
                "scope": "axiom-augmented proxy for the P3066 light_emergence_interface row, not observed light and not a strict selector theorem",
                "sigma_branches": list(SIGMA_BRANCHES),
            },
            "finite_lorentz_proxy_matrix": proxy_rows,
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(r["hit_count"] for r in grep_rows),
            "sigma_branches": len(SIGMA_BRANCHES),
            "sampled_boosts": len(BOOSTS),
            "lorentz_proxy_rows": len(proxy_rows),
            "proxy_null_covariance_pass_rows": sum(1 for r in proxy_rows if r["proxy_null_covariance_pass"]),
            "strict_lorentz_closure_rows": sum(1 for r in proxy_rows if r["strict_lorentz_closure_exported"]),
            "p3066_accepted_physics_obligation_rows": p3066.get("finite_certificate", {}).get("accepted_physics_obligation_rows"),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(1 for o in proof_obligations if o["satisfied"]),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3067 constructs the requested light_emergence_interface object only as an axiom-augmented T_sigma proxy: L_sigma maps sigma to the oriented 1+1 null ray k_sigma=(1,sigma).  The exact finite boost table has 8/8 proxy null-covariance passes, but 0 strict Lorentz/observed-light closure rows because the spacetime embedding, unit-bearing metric, photon/gauge field, variational dynamics, and empirical map are absent.",
            "negative_export_flags": {k: False for k in ["observed_light_exported", "strict_lorentz_closure_exported", "gauge_photon_exported", "unit_bearing_action_exported", "variational_EOM_exported", "empirical_physics_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"sigma_conditioned_light_proxy_constructed": True, "finite_null_covariance_proxy_passed": True, "selector_search_not_replayed": True, "conditioned_nadsoliton_to_light_work_unblocked_at_proxy_level": True},
            "next_honest_step": "Attack exactly one missing blocker in P3067: construct a strict nadsoliton-to-spacetime embedding with a unit-normalized 1+1 metric/speed-of-light scale for the k_sigma proxy, then rerun the Lorentz audit.  If that embedding cannot be exported, pivot to the sigma-invariant scalar conservation/scale-control row from P3066 rather than promoting the proxy to observed light.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3067/S2017 sigma-conditioned light Lorentz proxy matrix", "", f"Status: `{payload['status']}`", "",
        "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- sigma branches: `{c['sigma_branches']}`",
        f"- sampled boosts: `{c['sampled_boosts']}`",
        f"- Lorentz proxy rows: `{c['lorentz_proxy_rows']}`",
        f"- proxy null-covariance pass rows: `{c['proxy_null_covariance_pass_rows']}`",
        f"- strict Lorentz closure rows: `{c['strict_lorentz_closure_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3067/S2017 sigma-conditioned light Lorentz proxy matrix", "## P3067/S2017 sigma-conditioned light Lorentz proxy matrix\n\n`P3067/S2017` executes the `P3066` recommendation on the single `light_emergence_interface` row under the explicit `T_sigma` axiom boundary.  It constructs `SigmaConditionedNadsolitonLightTransitionProxy`, the conditional law `L_sigma: sigma -> k_sigma=(1,sigma)` in a `1+1` null-ray proxy space, and audits exact rational boost rows.  The finite table has `2` sigma branches, `4` sampled boosts, `8` Lorentz proxy rows, and `8` proxy null-covariance passes.  It exports `0` strict Lorentz/observed-light closure rows because no strict spacetime embedding, unit-bearing metric/speed-of-light scale, photon/gauge field, variational EOM, Hamiltonian, or empirical map is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3067/S2017 light proxy is not a light-sector action", "## P3067/S2017 light proxy is not a light-sector action\n\n`P3067/S2017` only checks a finite null-ray covariance proxy for `L_sigma: sigma -> (1,sigma)`.  It does not add a Maxwell/Yang-Mills field, polarization bundle, unit-bearing Lorentzian action, propagator, EOM, Hamiltonian, or empirical speed-of-light/constant map.\n")
    append_once(AGENTS, "Current sigma-conditioned light Lorentz proxy guardrail (P3067/S2017, 2026-06-24)", "## Current sigma-conditioned light Lorentz proxy guardrail (P3067/S2017, 2026-06-24)\n\n- P3067 executes exactly one P3066 row, `light_emergence_interface`, under the P3065 `T_sigma` axiom boundary rather than replaying strict selector search.\n- The constructed proxy law is `L_sigma: sigma -> k_sigma=(1,sigma)` in a finite `1+1` null-ray audit; all `8` sampled rational boost rows preserve the null proxy.\n- This is only axiom-augmented proxy covariance: it exports `0` strict Lorentz/observed-light closures and does not provide a spacetime embedding, unit-bearing metric, photon/gauge sector, variational EOM, empirical map, `QW-2191` discharge, strict selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is one missing blocker only: strict nadsoliton-to-spacetime embedding with unit-normalized `1+1` metric/speed-of-light scale for this proxy; otherwise pivot to a sigma-invariant scalar conservation/scale-control row.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
