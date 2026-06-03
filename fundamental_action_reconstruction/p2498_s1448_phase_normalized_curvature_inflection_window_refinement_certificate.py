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
    interval_row,
)
from p2494_s1444_phase_normalized_compression_curvature_multiprecision_stability_certificate import (
    curvature_value,
    mp_params,
)

GEN = ROOT / "generated"
OUT = GEN / "p2498_s1448_phase_normalized_curvature_inflection_window_refinement_certificate.json"
MD = GEN / "p2498_s1448_phase_normalized_curvature_inflection_window_refinement_certificate.md"

SOURCE_FILES = {
    "P2497_INFLECTION_WINDOW": GEN / "p2497_s1447_phase_normalized_curvature_inflection_window_interval_isolation_certificate.json",
}

mp.mp.dps = 120
REFINED_LEFT = mp.mpf("0.34961674")
REFINED_RIGHT = mp.mpf("0.34961675")
MAX_DEPTH = 32
LEDGER_EDGE_ROWS = 4


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
        "new_packet": "P2498|S1448|inflection window refinement|root-window contraction|nested inflection window|submicro inflection|curvature root refinement",
        "precursor_packets": "P2493|S1443|P2494|S1444|P2495|S1445|P2496|S1446|P2497|S1447|inflection-window interval isolation",
        "refinement_language": "adaptive interval-slab|slab exclusion|root-window|bisection|window contraction|endpoint sign|curvature zero",
        "guardrail_language": "legacy -> strict completion bridge|role-transfer audit|silent inheritance|K_legacy_ont|K_strict_gate|QW-2191",
        "closure_blockers": "role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer|directed-rounding",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def text(value: mp.mpf, digits: int = 60) -> str:
    return mp.nstr(value, digits)


def point_curvature_sign(d: mp.mpf) -> int:
    mp.mp.dps = 120
    value = curvature_value(d, mp_params())["x_second"]
    return 1 if value > 0 else -1 if value < 0 else 0


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
    return {"accepted_cells": accepted, "unresolved_cells": unresolved}


def summarize_region(raw: dict[str, Any], left: mp.mpf, right: mp.mpf, expected_sign: int) -> dict[str, Any]:
    accepted = raw["accepted_cells"]
    unresolved = raw["unresolved_cells"]
    widths = [mp.mpf(row["width_d"]) for row in accepted]
    ledger_digest = sha256_json({"accepted_cells": accepted, "unresolved_cells": unresolved})
    return {
        "region_left_d": text(left, 30),
        "region_right_d": text(right, 30),
        "expected_sign": expected_sign,
        "accepted_cell_count": len(accepted),
        "unresolved_cell_count": len(unresolved),
        "max_subdivision_depth": max((row["subdivision_depth"] for row in accepted + unresolved), default=0),
        "max_accepted_cell_width": text(max(widths), 30) if widths else None,
        "min_accepted_cell_width": text(min(widths), 30) if widths else None,
        "all_cells_exclude_zero_with_expected_sign": not unresolved and all(row["curvature_interval_sign"] == expected_sign for row in accepted),
        "ledger_cell_count": len(accepted) + len(unresolved),
        "ledger_digest_sha256": ledger_digest,
        "accepted_cell_head": accepted[:LEDGER_EDGE_ROWS],
        "accepted_cell_tail": accepted[-LEDGER_EDGE_ROWS:],
        "unresolved_cells": unresolved,
    }


def build_refinement_certificate(p2497: dict[str, Any]) -> dict[str, Any]:
    prior = p2497["inflection_window_interval_isolation_certificate"]
    prior_left = mp.mpf(prior["inflection_root_window_d"][0])
    prior_right = mp.mpf(prior["inflection_root_window_d"][1])
    left_raw = certify_region(prior_left, REFINED_LEFT, 1)
    right_raw = certify_region(REFINED_RIGHT, prior_right, -1)
    left_summary = summarize_region(left_raw, prior_left, REFINED_LEFT, 1)
    right_summary = summarize_region(right_raw, REFINED_RIGHT, prior_right, -1)
    root_estimate = mp.mpf(prior["point_root_estimate"])
    refined_width = REFINED_RIGHT - REFINED_LEFT
    prior_width = prior_right - prior_left
    return {
        "prior_p2497_window_d": [text(prior_left, 30), text(prior_right, 30)],
        "prior_p2497_window_width": text(prior_width, 30),
        "refined_inflection_window_d": [text(REFINED_LEFT, 30), text(REFINED_RIGHT, 30)],
        "refined_inflection_window_width": text(refined_width, 30),
        "window_width_contraction_factor": text(prior_width / refined_width, 30),
        "left_refinement_positive_slab_exclusion": left_summary,
        "right_refinement_negative_slab_exclusion": right_summary,
        "point_root_estimate": text(root_estimate, 80),
        "point_root_inside_refined_window": REFINED_LEFT < root_estimate < REFINED_RIGHT,
        "refined_left_endpoint_point_curvature_sign": point_curvature_sign(REFINED_LEFT),
        "refined_right_endpoint_point_curvature_sign": point_curvature_sign(REFINED_RIGHT),
        "p2497_outside_window_exclusion_inherited": p2497.get("outside_window_interval_zero_exclusion_certified") is True,
        "outside_refined_window_zero_exclusion_certified_on_audited_domain": p2497.get("outside_window_interval_zero_exclusion_certified") is True and left_summary["all_cells_exclude_zero_with_expected_sign"] and right_summary["all_cells_exclude_zero_with_expected_sign"],
        "single_refined_unresolved_inflection_window_on_audited_domain": p2497.get("outside_window_interval_zero_exclusion_certified") is True and left_summary["all_cells_exclude_zero_with_expected_sign"] and right_summary["all_cells_exclude_zero_with_expected_sign"] and REFINED_LEFT < root_estimate < REFINED_RIGHT,
        "max_subdivision_depth_allowed": MAX_DEPTH,
        "stored_full_cell_ledger": False,
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
    body = path.read_text(encoding="utf-8")
    if marker not in body:
        path.write_text(body.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2498/S1448 phase-normalized curvature inflection-window refinement certificate

`P2498/S1448` contracts the P2497 unresolved inflection window from `[0.3495, 0.3498]` to `[0.34961674, 0.34961675]`.  Inside the old P2497 window, adaptive `mpmath.iv` slab exclusions certify positive curvature on the left shoulder and negative curvature on the right shoulder, while the P2497 outside-window exclusion is inherited for the rest of the audited domain.  The point root estimate remains inside the refined window, giving a much narrower finite localization of the compression-curvature sign flip.

This remains a finite interval-backed refinement certificate, not a formal directed-rounding proof backend, global analytic inflection uniqueness theorem, or nonlinear compression-flow source theorem.  It exports no damping bridge atom, no strict compression source, no selector/source closure, no QW-2191 discharge, no role-transfer, no role-bearing `L_total`, no physical-value generation, and no ToE closure.
"""
    lag_section = """
## P2498/S1448 curvature inflection-window refinement guard

`P2498/S1448` refines the P2497 localized sign-change window to `[0.34961674, 0.34961675]` by interval-excluding zero on both shoulders of the old window.  This sharpens the nonaffine compression-flow obstruction without promoting it to a source theorem, bridge atom, role-transfer theorem, QW-2191 discharge, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2498/S1448 phase-normalized curvature inflection-window refinement certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2498/S1448 curvature inflection-window refinement guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2497 = theorem(sources["P2497_INFLECTION_WINDOW"], "phase_normalized_curvature_inflection_window_interval_isolation_certificate")
    cert = build_refinement_certificate(p2497)
    theorem_export = {
        "theorem_name": "P2498_T1_phase_normalized_curvature_inflection_window_refinement_certificate",
        "audited_chain": ["P2493/S1443", "P2494/S1444", "P2495/S1445", "P2496/S1446", "P2497/S1447"],
        "inflection_window_refinement_certificate": cert,
        "prior_p2497_window_d": cert["prior_p2497_window_d"],
        "refined_inflection_window_d": cert["refined_inflection_window_d"],
        "refined_inflection_window_width": cert["refined_inflection_window_width"],
        "window_width_contraction_factor": cert["window_width_contraction_factor"],
        "point_root_inside_refined_window": cert["point_root_inside_refined_window"],
        "left_refinement_positive_cell_count": cert["left_refinement_positive_slab_exclusion"]["accepted_cell_count"],
        "right_refinement_negative_cell_count": cert["right_refinement_negative_slab_exclusion"]["accepted_cell_count"],
        "left_refinement_zero_excluded": cert["left_refinement_positive_slab_exclusion"]["all_cells_exclude_zero_with_expected_sign"],
        "right_refinement_zero_excluded": cert["right_refinement_negative_slab_exclusion"]["all_cells_exclude_zero_with_expected_sign"],
        "p2497_outside_window_exclusion_inherited": cert["p2497_outside_window_exclusion_inherited"],
        "outside_refined_window_zero_exclusion_certified_on_audited_domain": cert["outside_refined_window_zero_exclusion_certified_on_audited_domain"],
        "single_refined_unresolved_inflection_window_on_audited_domain": cert["single_refined_unresolved_inflection_window_on_audited_domain"],
        "stored_full_cell_ledger": cert["stored_full_cell_ledger"],
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
            "P2498 is a finite mpmath.iv window-refinement audit, not a formal directed-rounding proof backend.",
            "A narrower inflection window does not derive the nonlinear compression-flow law from strict information geometry.",
            "No damping bridge atom, selector/source theorem, QW-2191 discharge, role-transfer theorem, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Replace the mpmath.iv refinement with a formal directed-rounding/analytic monotonicity proof, or derive the nonlinear compression-flow source that explains the localized inflection.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2497_exclusion_inherited": theorem_export["p2497_outside_window_exclusion_inherited"],
        "refined_window_is_tiny": mp.mpf(theorem_export["refined_inflection_window_width"]) <= mp.mpf("1e-8"),
        "contraction_factor_large": mp.mpf(theorem_export["window_width_contraction_factor"]) >= mp.mpf("30000"),
        "point_root_inside_refined_window": theorem_export["point_root_inside_refined_window"],
        "left_refinement_zero_excluded": theorem_export["left_refinement_zero_excluded"],
        "right_refinement_zero_excluded": theorem_export["right_refinement_zero_excluded"],
        "outside_refined_window_zero_excluded": theorem_export["outside_refined_window_zero_exclusion_certified_on_audited_domain"],
        "single_refined_window_recorded": theorem_export["single_refined_unresolved_inflection_window_on_audited_domain"],
        "full_cell_ledger_not_stored_but_hashed": theorem_export["stored_full_cell_ledger"] is False and bool(cert["left_refinement_positive_slab_exclusion"]["ledger_digest_sha256"]) and bool(cert["right_refinement_negative_slab_exclusion"]["ledger_digest_sha256"]),
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
        "packet_id": "P2498",
        "stage_id": "S1448",
        "status": "PHASE_NORMALIZED_CURVATURE_INFLECTION_WINDOW_REFINEMENT_CERTIFICATE_NO_FORMAL_DIRECTED_BACKEND_NO_GLOBAL_UNIQUENESS_NO_SOURCE_EXPORT_NO_BRIDGE_ATOM_NO_QW2191_NO_ROLE_TRANSFER_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "phase_normalized_curvature_inflection_window_refinement_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["phase_normalized_curvature_inflection_window_refinement_certificate"]["theorem_export"]
    c = t["inflection_window_refinement_certificate"]
    lines = [
        "# P2498/S1448 phase-normalized curvature inflection-window refinement certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Prior P2497 window: `{c['prior_p2497_window_d']}` with width `{c['prior_p2497_window_width']}`.",
        f"- Refined unresolved window: `{c['refined_inflection_window_d']}` with width `{c['refined_inflection_window_width']}`.",
        f"- Width contraction factor: `{c['window_width_contraction_factor']}`.",
        f"- Point root estimate: `{c['point_root_estimate']}`; inside refined window: `{c['point_root_inside_refined_window']}`.",
        f"- Left shoulder cells: `{t['left_refinement_positive_cell_count']}`; zero excluded: `{t['left_refinement_zero_excluded']}`.",
        f"- Right shoulder cells: `{t['right_refinement_negative_cell_count']}`; zero excluded: `{t['right_refinement_zero_excluded']}`.",
        f"- P2497 outside-window exclusion inherited: `{t['p2497_outside_window_exclusion_inherited']}`.",
        "",
        "## Ledger digests",
        "",
        f"- Left shoulder digest: `{c['left_refinement_positive_slab_exclusion']['ledger_digest_sha256']}`.",
        f"- Right shoulder digest: `{c['right_refinement_negative_slab_exclusion']['ledger_digest_sha256']}`.",
        "",
        "## Negative controls",
        "",
        "This packet does not export a formal directed-rounding backend, global inflection uniqueness theorem, analytic monotonicity theorem, curvature dynamic source, legacy-to-strict bridge atom, strict compression source, selector/source theorem, QW-2191 discharge, role-transfer license, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['phase_normalized_curvature_inflection_window_refinement_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["phase_normalized_curvature_inflection_window_refinement_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
