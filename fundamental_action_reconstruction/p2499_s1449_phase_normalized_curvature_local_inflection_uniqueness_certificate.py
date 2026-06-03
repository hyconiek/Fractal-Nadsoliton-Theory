#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

import mpmath as mp

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)
from p2494_s1444_phase_normalized_compression_curvature_multiprecision_stability_certificate import (
    branch_x,
    mp_params,
)
from p2495_s1445_phase_normalized_compression_curvature_interval_enclosure_certificate import (
    interval_bounds_text,
    interval_parameters,
    interval_sign,
    interval_width_text,
    iv_interval,
    legacy_derivative_interval,
    legacy_second_interval,
    strict_derivative_interval,
    strict_second_interval,
)
from p2497_s1447_phase_normalized_curvature_inflection_window_interval_isolation_certificate import (
    interval_row,
)

GEN = ROOT / "generated"
OUT = GEN / "p2499_s1449_phase_normalized_curvature_local_inflection_uniqueness_certificate.json"
MD = GEN / "p2499_s1449_phase_normalized_curvature_local_inflection_uniqueness_certificate.md"

SOURCE_FILES = {
    "P2498_REFINED_WINDOW": GEN / "p2498_s1448_phase_normalized_curvature_inflection_window_refinement_certificate.json",
}

mp.mp.dps = 120
IV_DPS = 90
POINT_DPS = 120


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:40]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2499|S1449|local inflection uniqueness|curvature monotonicity window|third derivative interval|inflection uniqueness window|curvature derivative interval",
        "precursor_packets": "P2497|S1447|P2498|S1448|inflection-window interval isolation|inflection window refinement|root-window contraction",
        "third_derivative_language": "x'''|third derivative|curvature derivative|monotonicity window|local uniqueness|strictly decreasing curvature",
        "guardrail_language": "legacy -> strict completion bridge|role-transfer audit|silent inheritance|K_legacy_ont|K_strict_gate|QW-2191",
        "closure_blockers": "role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer|directed-rounding",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 60) -> str:
    return mp.nstr(value, digits)


def legacy_third_interval(x: mp.iv.mpf, p: dict[str, mp.iv.mpf]) -> mp.iv.mpf:
    omega = p["legacy_omega"]
    phi = p["legacy_phi"]
    beta = p["legacy_beta_tors"]
    theta = omega * x + phi
    cos_phi = mp.iv.cos(phi)
    f0 = mp.iv.cos(theta) / cos_phi
    f1 = -omega * mp.iv.sin(theta) / cos_phi
    f2 = -(omega**2) * mp.iv.cos(theta) / cos_phi
    f3 = (omega**3) * mp.iv.sin(theta) / cos_phi
    denominator = 1 + beta * x
    g0 = 1 / denominator
    g1 = -beta / denominator**2
    g2 = 2 * beta**2 / denominator**3
    g3 = -6 * beta**3 / denominator**4
    return f3 * g0 + 3 * f2 * g1 + 3 * f1 * g2 + f0 * g3


def strict_third_interval(d: mp.iv.mpf, p: dict[str, mp.iv.mpf]) -> mp.iv.mpf:
    omega = p["strict_omega"]
    phi = p["strict_phi"]
    beta = p["strict_beta"]
    eta = p["strict_eta"]
    theta = omega * d + phi
    cos_phi = mp.iv.cos(phi)
    f0 = mp.iv.cos(theta) / cos_phi
    f1 = -omega * mp.iv.sin(theta) / cos_phi
    f2 = -(omega**2) * mp.iv.cos(theta) / cos_phi
    f3 = (omega**3) * mp.iv.sin(theta) / cos_phi
    u0 = 1 + beta * d**eta
    u1 = beta * eta * d ** (eta - 1)
    u2 = beta * eta * (eta - 1) * d ** (eta - 2)
    u3 = beta * eta * (eta - 1) * (eta - 2) * d ** (eta - 3)
    h0 = 1 / u0
    h1 = -u1 / u0**2
    h2 = -u2 / u0**2 + 2 * u1**2 / u0**3
    h3 = -u3 / u0**2 + 6 * u1 * u2 / u0**3 - 6 * u1**3 / u0**4
    return f3 * h0 + 3 * f2 * h1 + 3 * f1 * h2 + f0 * h3


def third_derivative_interval(left: mp.mpf, right: mp.mpf) -> dict[str, Any]:
    mp.iv.dps = IV_DPS
    mp.mp.dps = POINT_DPS
    p_point = mp_params()
    p_interval = interval_parameters()
    x_left = branch_x(left, p_point)
    x_right = branch_x(right, p_point)
    x_interval = iv_interval(min(x_left, x_right), max(x_left, x_right))
    d_interval = mp.iv.mpf([left, right])
    legacy_prime = legacy_derivative_interval(x_interval, p_interval)
    legacy_second = legacy_second_interval(x_interval, p_interval)
    legacy_third = legacy_third_interval(x_interval, p_interval)
    strict_prime = strict_derivative_interval(d_interval, p_interval)
    strict_second = strict_second_interval(d_interval, p_interval)
    strict_third = strict_third_interval(d_interval, p_interval)
    x_prime = strict_prime / legacy_prime
    x_second = (strict_second - legacy_second * x_prime**2) / legacy_prime
    x_third = (strict_third - legacy_third * x_prime**3 - 3 * legacy_second * x_prime * x_second) / legacy_prime
    return {
        "refined_window_d": [text(left, 30), text(right, 30)],
        "x_interval_from_endpoint_branch_solve": [text(min(x_left, x_right), 80), text(max(x_left, x_right), 80)],
        "legacy_prime_interval": interval_bounds_text(legacy_prime),
        "legacy_second_interval": interval_bounds_text(legacy_second),
        "legacy_third_interval": interval_bounds_text(legacy_third),
        "strict_prime_interval": interval_bounds_text(strict_prime),
        "strict_second_interval": interval_bounds_text(strict_second),
        "strict_third_interval": interval_bounds_text(strict_third),
        "x_prime_interval": interval_bounds_text(x_prime),
        "x_second_interval_over_window": interval_bounds_text(x_second),
        "x_third_interval_over_window": interval_bounds_text(x_third),
        "x_third_interval_width": interval_width_text(x_third),
        "x_third_interval_sign": interval_sign(x_third),
        "x_third_strictly_negative_on_refined_window": interval_sign(x_third) == -1,
    }


def build_uniqueness_certificate(p2498: dict[str, Any]) -> dict[str, Any]:
    refined_left, refined_right = [mp.mpf(value) for value in p2498["refined_inflection_window_d"]]
    left_endpoint = interval_row(refined_left, refined_left, 1)
    right_endpoint = interval_row(refined_right, refined_right, -1)
    third = third_derivative_interval(refined_left, refined_right)
    local_unique = (
        left_endpoint["curvature_interval_sign"] == 1
        and right_endpoint["curvature_interval_sign"] == -1
        and third["x_third_strictly_negative_on_refined_window"]
    )
    return {
        "refined_inflection_window_d": [text(refined_left, 30), text(refined_right, 30)],
        "refined_inflection_window_width": text(refined_right - refined_left, 30),
        "left_endpoint_curvature_interval": left_endpoint,
        "right_endpoint_curvature_interval": right_endpoint,
        "third_derivative_interval_certificate": third,
        "endpoint_intervals_have_opposite_strict_signs": left_endpoint["curvature_interval_sign"] == 1 and right_endpoint["curvature_interval_sign"] == -1,
        "x_second_strictly_decreasing_on_refined_window": third["x_third_strictly_negative_on_refined_window"],
        "finite_local_unique_curvature_zero_in_refined_window": local_unique,
        "p2498_outside_refined_window_zero_exclusion_inherited": p2498.get("outside_refined_window_zero_exclusion_certified_on_audited_domain") is True,
        "single_local_inflection_plus_outside_exclusion_on_audited_domain": local_unique and p2498.get("outside_refined_window_zero_exclusion_certified_on_audited_domain") is True,
        "formal_directed_rounding_backend_exported": False,
        "global_analytic_inflection_uniqueness_theorem_exported": False,
        "analytic_third_derivative_symbolic_proof_exported": False,
        "curvature_dynamic_source_exported": False,
        "legacy_to_strict_bridge_atom_exported": False,
        "strict_compression_dynamic_source_exported": False,
        "selector_source_theorem_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "toe_closure_exported": False,
    }


def append_once(path, marker: str, section: str) -> None:
    body = path.read_text(encoding="utf-8")
    if marker not in body:
        path.write_text(body.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2499/S1449 phase-normalized curvature local inflection-uniqueness certificate

`P2499/S1449` adds a third-derivative interval check on the P2498 refined window `[0.34961674, 0.34961675]`.  The endpoint `x''` intervals have opposite strict signs, and the computed interval for `x'''` is strictly negative across the refined window.  Within the finite audited interval backend, this upgrades the refined window from a mere sign-change bracket to a local uniqueness certificate for the curvature zero in that window, while inheriting the P2498 outside-window zero exclusion for the rest of the audited domain.

This is still a finite `mpmath.iv` certificate, not a formal directed-rounding backend, global analytic inflection theorem, or nonlinear compression-flow source theorem.  It exports no damping bridge atom, no strict compression source, no selector/source closure, no QW-2191 discharge, no role-transfer, no role-bearing `L_total`, no physical-value generation, and no ToE closure.
"""
    lag_section = """
## P2499/S1449 local inflection-uniqueness guard

`P2499/S1449` checks the third derivative on the P2498 refined window and certifies that `x''` is interval-decreasing there with opposite endpoint signs.  This turns the localized sign-change guard into a finite local uniqueness guard for the audited branch, without exporting a nonlinear compression-flow source, bridge atom, role-transfer theorem, QW-2191 discharge, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2499/S1449 phase-normalized curvature local inflection-uniqueness certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2499/S1449 local inflection-uniqueness guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2498 = theorem(sources["P2498_REFINED_WINDOW"], "phase_normalized_curvature_inflection_window_refinement_certificate")
    cert = build_uniqueness_certificate(p2498)
    third = cert["third_derivative_interval_certificate"]
    theorem_export = {
        "theorem_name": "P2499_T1_phase_normalized_curvature_local_inflection_uniqueness_certificate",
        "audited_chain": ["P2493/S1443", "P2494/S1444", "P2495/S1445", "P2496/S1446", "P2497/S1447", "P2498/S1448"],
        "local_inflection_uniqueness_certificate": cert,
        "refined_inflection_window_d": cert["refined_inflection_window_d"],
        "refined_inflection_window_width": cert["refined_inflection_window_width"],
        "endpoint_intervals_have_opposite_strict_signs": cert["endpoint_intervals_have_opposite_strict_signs"],
        "x_third_interval_over_refined_window": third["x_third_interval_over_window"],
        "x_third_interval_sign": third["x_third_interval_sign"],
        "x_second_strictly_decreasing_on_refined_window": cert["x_second_strictly_decreasing_on_refined_window"],
        "finite_local_unique_curvature_zero_in_refined_window": cert["finite_local_unique_curvature_zero_in_refined_window"],
        "p2498_outside_refined_window_zero_exclusion_inherited": cert["p2498_outside_refined_window_zero_exclusion_inherited"],
        "single_local_inflection_plus_outside_exclusion_on_audited_domain": cert["single_local_inflection_plus_outside_exclusion_on_audited_domain"],
        "formal_directed_rounding_backend_exported": False,
        "global_analytic_inflection_uniqueness_theorem_exported": False,
        "analytic_third_derivative_symbolic_proof_exported": False,
        "curvature_dynamic_source_exported": False,
        "legacy_to_strict_bridge_atom_exported": False,
        "strict_compression_dynamic_source_exported": False,
        "selector_source_theorem_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "P2499 is a finite mpmath.iv local monotonicity certificate, not a formal directed-rounding proof backend.",
            "Local uniqueness of the audited curvature zero does not derive the nonlinear compression-flow law or a damping bridge atom.",
            "No selector/source theorem, QW-2191 discharge, role-transfer theorem, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Promote the third-derivative interval check to a formal directed-rounding or symbolic proof, or derive the strict nonlinear compression-flow source behind the locally unique inflection.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2498_outside_exclusion_inherited": theorem_export["p2498_outside_refined_window_zero_exclusion_inherited"],
        "endpoint_intervals_opposite_strict_signs": theorem_export["endpoint_intervals_have_opposite_strict_signs"],
        "third_derivative_interval_negative": theorem_export["x_third_interval_sign"] == -1,
        "x_second_strictly_decreasing": theorem_export["x_second_strictly_decreasing_on_refined_window"],
        "finite_local_unique_zero": theorem_export["finite_local_unique_curvature_zero_in_refined_window"],
        "single_local_plus_outside_exclusion": theorem_export["single_local_inflection_plus_outside_exclusion_on_audited_domain"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "formal_directed_rounding_backend_exported",
            "global_analytic_inflection_uniqueness_theorem_exported",
            "analytic_third_derivative_symbolic_proof_exported",
            "curvature_dynamic_source_exported",
            "legacy_to_strict_bridge_atom_exported",
            "strict_compression_dynamic_source_exported",
            "selector_source_theorem_exported",
            "qw2191_discharged_by_this_certificate",
            "role_transfer_licensed_by_this_certificate",
            "toe_closure_exported",
        ]),
    }
    return {
        "packet_id": "P2499",
        "stage_id": "S1449",
        "status": "PHASE_NORMALIZED_CURVATURE_LOCAL_INFLECTION_UNIQUENESS_CERTIFICATE_NO_FORMAL_DIRECTED_BACKEND_NO_GLOBAL_ANALYTIC_UNIQUENESS_NO_SOURCE_EXPORT_NO_BRIDGE_ATOM_NO_QW2191_NO_ROLE_TRANSFER_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "phase_normalized_curvature_local_inflection_uniqueness_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["phase_normalized_curvature_local_inflection_uniqueness_certificate"]["theorem_export"]
    c = t["local_inflection_uniqueness_certificate"]
    lines = [
        "# P2499/S1449 phase-normalized curvature local inflection-uniqueness certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Refined window: `{c['refined_inflection_window_d']}` with width `{c['refined_inflection_window_width']}`.",
        f"- Endpoint intervals have opposite strict signs: `{c['endpoint_intervals_have_opposite_strict_signs']}`.",
        f"- `x'''` interval over refined window: `{t['x_third_interval_over_refined_window']}` with sign `{t['x_third_interval_sign']}`.",
        f"- `x''` strictly decreases on refined window: `{t['x_second_strictly_decreasing_on_refined_window']}`.",
        f"- Finite local unique curvature zero in refined window: `{t['finite_local_unique_curvature_zero_in_refined_window']}`.",
        f"- P2498 outside refined-window exclusion inherited: `{t['p2498_outside_refined_window_zero_exclusion_inherited']}`.",
        "",
        "## Negative controls",
        "",
        "This packet does not export a formal directed-rounding backend, global analytic inflection uniqueness theorem, analytic symbolic third-derivative proof, curvature dynamic source, legacy-to-strict bridge atom, strict compression source, selector/source theorem, QW-2191 discharge, role-transfer license, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['phase_normalized_curvature_local_inflection_uniqueness_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["phase_normalized_curvature_local_inflection_uniqueness_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
