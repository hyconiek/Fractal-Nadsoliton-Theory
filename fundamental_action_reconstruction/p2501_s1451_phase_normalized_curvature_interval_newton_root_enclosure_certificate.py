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
from p2497_s1447_phase_normalized_curvature_inflection_window_interval_isolation_certificate import (
    curvature_interval_for_slab,
    interval_row,
)
from p2499_s1449_phase_normalized_curvature_local_inflection_uniqueness_certificate import (
    third_derivative_interval,
)

GEN = ROOT / "generated"
OUT = GEN / "p2501_s1451_phase_normalized_curvature_interval_newton_root_enclosure_certificate.json"
MD = GEN / "p2501_s1451_phase_normalized_curvature_interval_newton_root_enclosure_certificate.md"

SOURCE_FILES = {
    "P2499_LOCAL_UNIQUENESS": GEN / "p2499_s1449_phase_normalized_curvature_local_inflection_uniqueness_certificate.json",
    "P2500_SYMBOLIC_PROVENANCE": GEN / "p2500_s1450_phase_normalized_third_derivative_symbolic_identity_audit.json",
}

mp.mp.dps = 120
mp.iv.dps = 90
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
        "new_packet": "P2501|S1451|interval Newton root enclosure|Newton contraction|contractive root enclosure|curvature zero enclosure|inflection root enclosure",
        "precursor_packets": "P2497|S1447|P2498|S1448|P2499|S1449|P2500|S1450|inflection window|local inflection uniqueness|third-derivative symbolic identity",
        "interval_newton_language": "interval Newton|Newton enclosure|Krawczyk|contractive root enclosure|self-mapping interval|root residual|root enclosure",
        "guardrail_language": "legacy -> strict completion bridge|role-transfer audit|silent inheritance|K_legacy_ont|K_strict_gate|QW-2191",
        "closure_blockers": "role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer|directed-rounding",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 60) -> str:
    return mp.nstr(value, digits)


def interval_from_text(bounds: list[str]) -> mp.iv.mpf:
    return mp.iv.mpf([mp.mpf(bounds[0]), mp.mpf(bounds[1])])


def interval_width(bounds: list[str]) -> mp.mpf:
    return mp.mpf(bounds[1]) - mp.mpf(bounds[0])


def interval_newton_step(left: mp.mpf, right: mp.mpf) -> dict[str, Any]:
    mp.mp.dps = POINT_DPS
    mp.iv.dps = IV_DPS
    midpoint = (left + right) / 2
    f_mid = curvature_interval_for_slab(midpoint, midpoint)
    derivative_cert = third_derivative_interval(left, right)
    derivative_interval = interval_from_text(derivative_cert["x_third_interval_over_window"])
    midpoint_interval = mp.iv.mpf([midpoint, midpoint])
    newton_image = midpoint_interval - f_mid / derivative_interval
    newton_left = mp.mpf(newton_image.a)
    newton_right = mp.mpf(newton_image.b)
    original_width = right - left
    newton_width = newton_right - newton_left
    left_endpoint = interval_row(newton_left, newton_left, 1)
    right_endpoint = interval_row(newton_right, newton_right, -1)
    return {
        "starting_window_d": [text(left, 40), text(right, 40)],
        "starting_window_width": text(original_width, 40),
        "midpoint_d": text(midpoint, 80),
        "curvature_interval_at_midpoint": [text(mp.mpf(f_mid.a), 90), text(mp.mpf(f_mid.b), 90)],
        "curvature_midpoint_interval_sign": 0 if mp.mpf(f_mid.a) <= 0 <= mp.mpf(f_mid.b) else (1 if mp.mpf(f_mid.a) > 0 else -1),
        "third_derivative_interval_over_starting_window": derivative_cert["x_third_interval_over_window"],
        "third_derivative_interval_sign": derivative_cert["x_third_interval_sign"],
        "third_derivative_strictly_negative": derivative_cert["x_third_strictly_negative_on_refined_window"],
        "interval_newton_image_d": [text(newton_left, 90), text(newton_right, 90)],
        "interval_newton_image_width": text(newton_width, 90),
        "contraction_factor": text(original_width / newton_width, 50),
        "newton_image_subset_of_starting_window": left <= newton_left and newton_right <= right,
        "newton_image_strict_subset_of_starting_window": left < newton_left and newton_right < right,
        "left_endpoint_curvature_interval": left_endpoint,
        "right_endpoint_curvature_interval": right_endpoint,
        "contracted_endpoint_intervals_have_opposite_strict_signs": left_endpoint["curvature_interval_sign"] == 1 and right_endpoint["curvature_interval_sign"] == -1,
        "contracted_enclosure_contains_a_curvature_zero_by_endpoint_signs": left_endpoint["curvature_interval_sign"] == 1 and right_endpoint["curvature_interval_sign"] == -1,
        "contracted_enclosure_contains_unique_curvature_zero_by_p2499_monotonicity": False,
    }


def build_newton_certificate(p2499: dict[str, Any], p2500: dict[str, Any]) -> dict[str, Any]:
    left, right = [mp.mpf(value) for value in p2499["refined_inflection_window_d"]]
    step = interval_newton_step(left, right)
    unique = (
        step["newton_image_subset_of_starting_window"]
        and step["third_derivative_strictly_negative"]
        and step["contracted_endpoint_intervals_have_opposite_strict_signs"]
        and p2499.get("finite_local_unique_curvature_zero_in_refined_window") is True
    )
    step["contracted_enclosure_contains_unique_curvature_zero_by_p2499_monotonicity"] = unique
    return {
        "interval_backend": "mpmath.iv",
        "interval_dps": IV_DPS,
        "point_dps": POINT_DPS,
        "interval_newton_formula": "N(m,I)=m - x_second(m)/x_third(I)",
        "interval_newton_step": step,
        "p2499_local_uniqueness_inherited": p2499.get("finite_local_unique_curvature_zero_in_refined_window") is True,
        "p2499_outside_refined_window_zero_exclusion_inherited": p2499.get("single_local_inflection_plus_outside_exclusion_on_audited_domain") is True,
        "p2500_symbolic_formula_provenance_inherited": p2500.get("all_symbolic_residuals_zero") is True,
        "contractive_interval_newton_enclosure_certified": step["newton_image_subset_of_starting_window"] and step["third_derivative_strictly_negative"],
        "unique_inflection_root_in_contracted_enclosure_on_audited_branch": unique,
        "single_inflection_root_on_audited_domain_inherited_and_narrowed": unique and p2499.get("single_local_inflection_plus_outside_exclusion_on_audited_domain") is True,
        "formal_directed_rounding_backend_exported": False,
        "global_analytic_inflection_uniqueness_theorem_exported": False,
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
## P2501/S1451 phase-normalized curvature interval-Newton root enclosure certificate

`P2501/S1451` applies a finite interval-Newton contraction to the P2498/P2499 refined inflection window.  With `f(d)=x''(d)`, midpoint `m=0.349616745`, and the P2499/P2500-audited negative interval for `f'(d)=x'''(d)`, the image `m - f(m)/f'([0.34961674,0.34961675])` is a strict subset of the starting window, approximately `[0.3496167445840840099, 0.3496167445840840974]`.  Endpoint curvature signs on the contracted interval remain opposite, and P2499 monotonicity supplies local uniqueness inside the contracted enclosure.

This is a finite `mpmath.iv` interval-Newton enclosure and formula-audited numerical contraction, not a formal directed-rounding backend, global analytic inflection theorem, or nonlinear compression-flow source theorem.  It exports no damping bridge atom, no strict compression source, no selector/source closure, no QW-2191 discharge, no role-transfer, no role-bearing `L_total`, no physical-value generation, and no ToE closure.
"""
    lag_section = """
## P2501/S1451 interval-Newton inflection-root enclosure guard

`P2501/S1451` contracts the P2498/P2499 local inflection bracket by an interval-Newton step using the formula-audited `x'''` denominator.  The contracted enclosure is a strict subset of the refined window and inherits local uniqueness from P2499, but it still exports no nonlinear compression-flow source, bridge atom, role-transfer theorem, QW-2191 discharge, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2501/S1451 phase-normalized curvature interval-Newton root enclosure certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2501/S1451 interval-Newton inflection-root enclosure guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2499 = theorem(sources["P2499_LOCAL_UNIQUENESS"], "phase_normalized_curvature_local_inflection_uniqueness_certificate")
    p2500 = theorem(sources["P2500_SYMBOLIC_PROVENANCE"], "phase_normalized_third_derivative_symbolic_identity_audit")
    cert = build_newton_certificate(p2499, p2500)
    step = cert["interval_newton_step"]
    theorem_export = {
        "theorem_name": "P2501_T1_phase_normalized_curvature_interval_newton_root_enclosure_certificate",
        "audited_chain": ["P2497/S1447", "P2498/S1448", "P2499/S1449", "P2500/S1450"],
        "interval_newton_root_enclosure_certificate": cert,
        "starting_window_d": step["starting_window_d"],
        "starting_window_width": step["starting_window_width"],
        "contracted_root_enclosure_d": step["interval_newton_image_d"],
        "contracted_root_enclosure_width": step["interval_newton_image_width"],
        "contraction_factor": step["contraction_factor"],
        "newton_image_subset_of_starting_window": step["newton_image_subset_of_starting_window"],
        "newton_image_strict_subset_of_starting_window": step["newton_image_strict_subset_of_starting_window"],
        "third_derivative_interval_sign": step["third_derivative_interval_sign"],
        "third_derivative_strictly_negative": step["third_derivative_strictly_negative"],
        "contracted_endpoint_intervals_have_opposite_strict_signs": step["contracted_endpoint_intervals_have_opposite_strict_signs"],
        "contractive_interval_newton_enclosure_certified": cert["contractive_interval_newton_enclosure_certified"],
        "unique_inflection_root_in_contracted_enclosure_on_audited_branch": cert["unique_inflection_root_in_contracted_enclosure_on_audited_branch"],
        "single_inflection_root_on_audited_domain_inherited_and_narrowed": cert["single_inflection_root_on_audited_domain_inherited_and_narrowed"],
        "p2499_local_uniqueness_inherited": cert["p2499_local_uniqueness_inherited"],
        "p2499_outside_refined_window_zero_exclusion_inherited": cert["p2499_outside_refined_window_zero_exclusion_inherited"],
        "p2500_symbolic_formula_provenance_inherited": cert["p2500_symbolic_formula_provenance_inherited"],
        "formal_directed_rounding_backend_exported": False,
        "global_analytic_inflection_uniqueness_theorem_exported": False,
        "curvature_dynamic_source_exported": False,
        "legacy_to_strict_bridge_atom_exported": False,
        "strict_compression_dynamic_source_exported": False,
        "selector_source_theorem_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "P2501 contracts a local root enclosure; it is not a formal directed-rounding proof backend.",
            "The contracted inflection root is still a compression-geometry diagnostic, not the nonlinear compression-flow source or a damping bridge atom.",
            "No selector/source theorem, QW-2191 discharge, role-transfer theorem, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Either replace mpmath.iv with a formal directed-rounding backend, or use the now sharply enclosed inflection geometry to propose and audit a real nonlinear compression-flow source/bridge atom.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2499_local_uniqueness_inherited": theorem_export["p2499_local_uniqueness_inherited"],
        "p2500_symbolic_formula_provenance_inherited": theorem_export["p2500_symbolic_formula_provenance_inherited"],
        "third_derivative_strictly_negative": theorem_export["third_derivative_strictly_negative"],
        "newton_image_subset_of_starting_window": theorem_export["newton_image_subset_of_starting_window"],
        "newton_image_strict_subset_of_starting_window": theorem_export["newton_image_strict_subset_of_starting_window"],
        "contracted_endpoint_signs_opposite": theorem_export["contracted_endpoint_intervals_have_opposite_strict_signs"],
        "contracted_width_smaller_than_starting_width": mp.mpf(theorem_export["contracted_root_enclosure_width"]) < mp.mpf(theorem_export["starting_window_width"]),
        "contraction_factor_large": mp.mpf(theorem_export["contraction_factor"]) > mp.mpf("10000000"),
        "single_inflection_root_inherited_and_narrowed": theorem_export["single_inflection_root_on_audited_domain_inherited_and_narrowed"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "formal_directed_rounding_backend_exported",
            "global_analytic_inflection_uniqueness_theorem_exported",
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
        "packet_id": "P2501",
        "stage_id": "S1451",
        "status": "PHASE_NORMALIZED_CURVATURE_INTERVAL_NEWTON_ROOT_ENCLOSURE_CERTIFICATE_NO_FORMAL_DIRECTED_BACKEND_NO_GLOBAL_ANALYTIC_UNIQUENESS_NO_SOURCE_EXPORT_NO_BRIDGE_ATOM_NO_QW2191_NO_ROLE_TRANSFER_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "phase_normalized_curvature_interval_newton_root_enclosure_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["phase_normalized_curvature_interval_newton_root_enclosure_certificate"]["theorem_export"]
    lines = [
        "# P2501/S1451 phase-normalized curvature interval-Newton root enclosure certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Starting P2498/P2499 window: `{t['starting_window_d']}` with width `{t['starting_window_width']}`.",
        f"- Contracted interval-Newton root enclosure: `{t['contracted_root_enclosure_d']}`.",
        f"- Contracted enclosure width: `{t['contracted_root_enclosure_width']}`.",
        f"- Contraction factor: `{t['contraction_factor']}`.",
        f"- Newton image subset of starting window: `{t['newton_image_subset_of_starting_window']}`.",
        f"- Third derivative strictly negative on starting window: `{t['third_derivative_strictly_negative']}`.",
        f"- Contracted endpoint curvature signs remain opposite: `{t['contracted_endpoint_intervals_have_opposite_strict_signs']}`.",
        f"- Unique inflection root in contracted enclosure on audited branch: `{t['unique_inflection_root_in_contracted_enclosure_on_audited_branch']}`.",
        "",
        "## Negative controls",
        "",
        "This packet does not export a formal directed-rounding backend, global analytic inflection uniqueness theorem, curvature dynamic source, legacy-to-strict bridge atom, strict compression source, selector/source theorem, QW-2191 discharge, role-transfer license, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['phase_normalized_curvature_interval_newton_root_enclosure_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["phase_normalized_curvature_interval_newton_root_enclosure_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
