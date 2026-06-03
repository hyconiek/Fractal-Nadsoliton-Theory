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
    curvature_value,
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

GEN = ROOT / "generated"
OUT = GEN / "p2497_s1447_phase_normalized_curvature_inflection_window_interval_isolation_certificate.json"
MD = GEN / "p2497_s1447_phase_normalized_curvature_inflection_window_interval_isolation_certificate.md"

SOURCE_FILES = {
    "P2496_INVERSE_BRANCH": GEN / "p2496_s1446_phase_normalized_inverse_branch_interval_existence_uniqueness_certificate.json",
}

mp.mp.dps = 120
IV_DPS = 90
POINT_DPS = 100
AUDITED_LEFT = mp.mpf("0.0001")
AUDITED_RIGHT = mp.mpf("11")
ROOT_WINDOW_LEFT = mp.mpf("0.3495")
ROOT_WINDOW_RIGHT = mp.mpf("0.3498")
MAX_DEPTH = 22


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
        "new_packet": "P2497|S1447|inflection window interval isolation|curvature inflection window|curvature root-window|interval slab exclusion|single unresolved inflection window",
        "precursor_packets": "P2493|S1443|P2494|S1444|P2495|S1445|P2496|S1446|compression curvature|inverse branch",
        "root_isolation_language": "inflection root|root isolation|interval Newton|sign-changing bracket|slab exclusion|curvature zero|monotonicity window",
        "guardrail_language": "legacy -> strict completion bridge|role-transfer audit|silent inheritance|K_legacy_ont|K_strict_gate|QW-2191",
        "closure_blockers": "role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer|root-window theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 60) -> str:
    return mp.nstr(value, digits)


def curvature_interval_for_slab(left: mp.mpf, right: mp.mpf) -> mp.iv.mpf:
    mp.iv.dps = IV_DPS
    mp.mp.dps = POINT_DPS
    p = mp_params()
    ip = interval_parameters()
    x_left = branch_x(left, p)
    x_right = branch_x(right, p)
    x_interval = iv_interval(min(x_left, x_right), max(x_left, x_right))
    d_interval = mp.iv.mpf([left, right])
    legacy_prime = legacy_derivative_interval(x_interval, ip)
    strict_prime = strict_derivative_interval(d_interval, ip)
    x_prime = strict_prime / legacy_prime
    return (strict_second_interval(d_interval, ip) - legacy_second_interval(x_interval, ip) * x_prime**2) / legacy_prime


def interval_row(left: mp.mpf, right: mp.mpf, expected_sign: int) -> dict[str, Any]:
    value = curvature_interval_for_slab(left, right)
    sign = interval_sign(value)
    return {
        "left_d": text(left, 30),
        "right_d": text(right, 30),
        "width_d": text(right - left, 30),
        "expected_sign": expected_sign,
        "curvature_interval": interval_bounds_text(value),
        "curvature_interval_width": interval_width_text(value),
        "curvature_interval_sign": sign,
        "strictly_excludes_zero_with_expected_sign": sign == expected_sign,
    }


def certify_region(left: mp.mpf, right: mp.mpf, expected_sign: int) -> dict[str, Any]:
    accepted: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    stack = [(left, right, 0)]
    while stack:
        a, b, depth = stack.pop()
        row = interval_row(a, b, expected_sign)
        if row["strictly_excludes_zero_with_expected_sign"]:
            row["subdivision_depth"] = depth
            accepted.append(row)
            continue
        if depth >= MAX_DEPTH:
            row["subdivision_depth"] = depth
            unresolved.append(row)
            continue
        mid = (a + b) / 2
        stack.append((mid, b, depth + 1))
        stack.append((a, mid, depth + 1))
    accepted.sort(key=lambda row: mp.mpf(row["left_d"]))
    unresolved.sort(key=lambda row: mp.mpf(row["left_d"]))
    widths = [mp.mpf(row["width_d"]) for row in accepted]
    return {
        "region_left_d": text(left, 30),
        "region_right_d": text(right, 30),
        "expected_sign": expected_sign,
        "accepted_cell_count": len(accepted),
        "unresolved_cell_count": len(unresolved),
        "max_subdivision_depth": max((row["subdivision_depth"] for row in accepted + unresolved), default=0),
        "max_accepted_cell_width": text(max(widths), 30) if widths else None,
        "min_accepted_cell_width": text(min(widths), 30) if widths else None,
        "all_cells_exclude_zero_with_expected_sign": not unresolved and all(row["strictly_excludes_zero_with_expected_sign"] for row in accepted),
        "accepted_cells": accepted,
        "unresolved_cells": unresolved,
    }


def endpoint_interval(d: mp.mpf) -> dict[str, Any]:
    row = interval_row(d, d, 1)
    row["point_curvature_sign"] = interval_sign(curvature_interval_for_slab(d, d))
    return row


def build_inflection_window_certificate(p2496: dict[str, Any]) -> dict[str, Any]:
    mp.mp.dps = POINT_DPS
    p = mp_params()
    root = mp.findroot(lambda z: curvature_value(z, p)["x_second"], (mp.mpf("0.3"), mp.mpf("0.4")))
    root_curvature_abs = abs(curvature_value(root, p)["x_second"])
    left_region = certify_region(AUDITED_LEFT, ROOT_WINDOW_LEFT, 1)
    right_region = certify_region(ROOT_WINDOW_RIGHT, AUDITED_RIGHT, -1)
    left_endpoint = endpoint_interval(ROOT_WINDOW_LEFT)
    right_endpoint = endpoint_interval(ROOT_WINDOW_RIGHT)
    left_endpoint["expected_sign"] = 1
    left_endpoint["strictly_excludes_zero_with_expected_sign"] = left_endpoint["curvature_interval_sign"] == 1
    right_endpoint["expected_sign"] = -1
    right_endpoint["strictly_excludes_zero_with_expected_sign"] = right_endpoint["curvature_interval_sign"] == -1
    window_width = ROOT_WINDOW_RIGHT - ROOT_WINDOW_LEFT
    return {
        "audited_domain_d": [text(AUDITED_LEFT, 30), text(AUDITED_RIGHT, 30)],
        "inflection_root_window_d": [text(ROOT_WINDOW_LEFT, 30), text(ROOT_WINDOW_RIGHT, 30)],
        "inflection_root_window_width": text(window_width, 30),
        "interval_backend": "mpmath.iv",
        "interval_dps": IV_DPS,
        "point_solve_dps": POINT_DPS,
        "max_subdivision_depth_allowed": MAX_DEPTH,
        "left_positive_slab_exclusion": left_region,
        "right_negative_slab_exclusion": right_region,
        "window_left_endpoint_interval": left_endpoint,
        "window_right_endpoint_interval": right_endpoint,
        "point_root_estimate": text(root, 80),
        "point_root_inside_window": ROOT_WINDOW_LEFT < root < ROOT_WINDOW_RIGHT,
        "point_root_curvature_abs": text(root_curvature_abs, 40),
        "endpoint_intervals_have_opposite_strict_signs": left_endpoint["curvature_interval_sign"] == 1 and right_endpoint["curvature_interval_sign"] == -1,
        "outside_window_interval_zero_exclusion_certified": left_region["all_cells_exclude_zero_with_expected_sign"] and right_region["all_cells_exclude_zero_with_expected_sign"],
        "single_unresolved_inflection_window_on_audited_domain": left_region["all_cells_exclude_zero_with_expected_sign"] and right_region["all_cells_exclude_zero_with_expected_sign"] and ROOT_WINDOW_LEFT < root < ROOT_WINDOW_RIGHT,
        "p2496_inverse_branch_inherited": p2496.get("all_sample_targets_inside_branch_range") is True and p2496.get("legacy_derivative_interval_strictly_negative") is True,
        "formal_directed_rounding_backend_exported": False,
        "global_inflection_uniqueness_theorem_exported": False,
        "analytic_monotonicity_theorem_exported": False,
        "curvature_dynamic_source_exported": False,
        "legacy_to_strict_bridge_atom_exported": False,
        "strict_compression_dynamic_source_exported": False,
        "selector_source_theorem_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "toe_closure_exported": False,
    }


def append_once(path, marker: str, section: str) -> None:
    text0 = path.read_text(encoding="utf-8")
    if marker not in text0:
        path.write_text(text0.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2497/S1447 phase-normalized curvature inflection-window interval isolation certificate

`P2497/S1447` upgrades the P2493--P2496 curvature audit from isolated sample signs to an adaptive interval-slab exclusion on the audited domain `d in [0.0001, 11]`.  The script encloses `x''(d)` on the positive side `[0.0001, 0.3495]` and the negative side `[0.3498, 11]` using the phase-normalized inverse branch and `mpmath.iv` derivative expressions.  Every accepted slab excludes zero with the expected sign, while the point root estimate lies inside the remaining window `[0.3495, 0.3498]`; the endpoint intervals of that window have opposite strict signs.

This is a finite adaptive interval-slab isolation certificate for the audited branch, not a formal directed-rounding proof backend, global analytic inflection uniqueness theorem, or nonlinear compression-flow source theorem.  It exports no damping bridge atom, no strict compression source, no selector/source closure, no QW-2191 discharge, no role-transfer, no role-bearing `L_total`, no physical-value generation, and no ToE closure.
"""
    lag_section = """
## P2497/S1447 curvature inflection-window isolation guard

`P2497/S1447` adds an adaptive interval-slab guard behind the nonaffine compression conclusion: outside the narrow audited window `[0.3495, 0.3498]`, the phase-normalized branch curvature is interval-excluded from zero with the expected sign.  This makes the observed curvature sign flip a localized compression-flow obstruction rather than a sparse-sample artifact, but it still does not export a nonlinear compression-flow source, bridge atom, role-transfer theorem, QW-2191 discharge, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2497/S1447 phase-normalized curvature inflection-window interval isolation certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2497/S1447 curvature inflection-window isolation guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2496 = theorem(sources["P2496_INVERSE_BRANCH"], "phase_normalized_inverse_branch_interval_existence_uniqueness_certificate")
    cert = build_inflection_window_certificate(p2496)
    theorem_export = {
        "theorem_name": "P2497_T1_phase_normalized_curvature_inflection_window_interval_isolation_certificate",
        "audited_chain": ["P2493/S1443", "P2494/S1444", "P2495/S1445", "P2496/S1446"],
        "inflection_window_interval_isolation_certificate": cert,
        "interval_backend": cert["interval_backend"],
        "interval_dps": cert["interval_dps"],
        "point_solve_dps": cert["point_solve_dps"],
        "inflection_root_window_d": cert["inflection_root_window_d"],
        "inflection_root_window_width": cert["inflection_root_window_width"],
        "point_root_estimate": cert["point_root_estimate"],
        "point_root_inside_window": cert["point_root_inside_window"],
        "left_positive_slab_cell_count": cert["left_positive_slab_exclusion"]["accepted_cell_count"],
        "right_negative_slab_cell_count": cert["right_negative_slab_exclusion"]["accepted_cell_count"],
        "outside_window_interval_zero_exclusion_certified": cert["outside_window_interval_zero_exclusion_certified"],
        "endpoint_intervals_have_opposite_strict_signs": cert["endpoint_intervals_have_opposite_strict_signs"],
        "single_unresolved_inflection_window_on_audited_domain": cert["single_unresolved_inflection_window_on_audited_domain"],
        "p2496_inverse_branch_inherited": cert["p2496_inverse_branch_inherited"],
        "formal_directed_rounding_backend_exported": False,
        "global_inflection_uniqueness_theorem_exported": False,
        "analytic_monotonicity_theorem_exported": False,
        "curvature_dynamic_source_exported": False,
        "legacy_to_strict_bridge_atom_exported": False,
        "strict_compression_dynamic_source_exported": False,
        "selector_source_theorem_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "P2497 is an adaptive mpmath.iv slab-exclusion audit, not a formal directed-rounding proof backend.",
            "The remaining inflection window is localized on the audited branch, but no global analytic uniqueness theorem is exported.",
            "No nonlinear compression-flow source, damping bridge atom, selector/source theorem, QW-2191 discharge, role-transfer theorem, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Either replace the mpmath.iv slab audit with a formal directed-rounding backend/global monotonicity proof, or derive the nonlinear compression-flow source/bridge atom that explains the isolated inflection window.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2496_inverse_branch_inherited": theorem_export["p2496_inverse_branch_inherited"],
        "root_window_is_narrow": mp.mpf(theorem_export["inflection_root_window_width"]) <= mp.mpf("0.0003"),
        "point_root_inside_window": theorem_export["point_root_inside_window"],
        "endpoint_intervals_opposite_strict_signs": theorem_export["endpoint_intervals_have_opposite_strict_signs"],
        "outside_window_zero_excluded": theorem_export["outside_window_interval_zero_exclusion_certified"],
        "left_region_has_cells": theorem_export["left_positive_slab_cell_count"] > 0,
        "right_region_has_cells": theorem_export["right_negative_slab_cell_count"] > 0,
        "single_unresolved_window_recorded": theorem_export["single_unresolved_inflection_window_on_audited_domain"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "formal_directed_rounding_backend_exported",
            "global_inflection_uniqueness_theorem_exported",
            "analytic_monotonicity_theorem_exported",
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
        "packet_id": "P2497",
        "stage_id": "S1447",
        "status": "PHASE_NORMALIZED_CURVATURE_INFLECTION_WINDOW_INTERVAL_ISOLATION_CERTIFICATE_NO_FORMAL_DIRECTED_BACKEND_NO_GLOBAL_UNIQUENESS_NO_SOURCE_EXPORT_NO_BRIDGE_ATOM_NO_QW2191_NO_ROLE_TRANSFER_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "phase_normalized_curvature_inflection_window_interval_isolation_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["phase_normalized_curvature_inflection_window_interval_isolation_certificate"]["theorem_export"]
    c = t["inflection_window_interval_isolation_certificate"]
    lines = [
        "# P2497/S1447 phase-normalized curvature inflection-window interval isolation certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Audited domain: `{c['audited_domain_d']}`.",
        f"- Isolated unresolved inflection window: `{c['inflection_root_window_d']}` with width `{c['inflection_root_window_width']}`.",
        f"- Point root estimate: `{c['point_root_estimate']}`; inside window: `{c['point_root_inside_window']}`.",
        f"- Left positive slab cells: `{t['left_positive_slab_cell_count']}`; all exclude zero: `{c['left_positive_slab_exclusion']['all_cells_exclude_zero_with_expected_sign']}`.",
        f"- Right negative slab cells: `{t['right_negative_slab_cell_count']}`; all exclude zero: `{c['right_negative_slab_exclusion']['all_cells_exclude_zero_with_expected_sign']}`.",
        f"- Endpoint intervals have opposite strict signs: `{t['endpoint_intervals_have_opposite_strict_signs']}`.",
        f"- P2496 inverse-branch guard inherited: `{t['p2496_inverse_branch_inherited']}`.",
        "",
        "## Negative controls",
        "",
        "This packet does not export a formal directed-rounding backend, global inflection uniqueness theorem, analytic monotonicity theorem, curvature dynamic source, legacy-to-strict bridge atom, strict compression source, selector/source theorem, QW-2191 discharge, role-transfer license, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['phase_normalized_curvature_inflection_window_interval_isolation_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["phase_normalized_curvature_inflection_window_interval_isolation_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
