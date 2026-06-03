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
    SAMPLE_D_STRINGS,
    Z12_D_STRINGS,
)

GEN = ROOT / "generated"
OUT = GEN / "p2496_s1446_phase_normalized_inverse_branch_interval_existence_uniqueness_certificate.json"
MD = GEN / "p2496_s1446_phase_normalized_inverse_branch_interval_existence_uniqueness_certificate.md"

SOURCE_FILES = {
    "P2495_INTERVAL_ENCLOSURE": GEN / "p2495_s1445_phase_normalized_compression_curvature_interval_enclosure_certificate.json",
}

IV_DPS = 100
POINT_DPS = 100
BRANCH_LEFT = mp.mpf("0")
BRANCH_RIGHT = mp.mpf("2")


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
        "new_packet": "P2496|S1446|inverse branch existence|branch uniqueness|phase-normalized inverse branch|legacy norm monotonic|inverse-branch existence|branch bracket",
        "precursor_packets": "P2493|S1443|P2494|S1444|P2495|S1445|curvature interval enclosure|nonaffine curvature",
        "branch_language": "L_norm\\(x\\(d\\)\\)|S_strict_norm|monotone inverse branch|branch bracket|legacy derivative|strict normalized output",
        "guardrail_language": "legacy -> strict completion bridge|role-transfer audit|silent inheritance|K_legacy_ont|K_strict_gate|QW-2191",
        "closure_blockers": "role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer|root-window theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def mp_params() -> dict[str, mp.mpf]:
    return {
        "legacy_omega": mp.pi / 4,
        "legacy_phi": mp.pi / 6,
        "legacy_beta_tors": mp.mpf("0.01"),
        "strict_omega": mp.mpf("0.18575"),
        "strict_phi": mp.mpf("0.16250"),
        "strict_beta": mp.mpf("1"),
        "strict_eta": mp.mpf("1.8"),
    }


def iv_mpf(value: mp.mpf | str) -> mp.iv.mpf:
    return mp.iv.mpf([value, value])


def interval_bounds_text(value: mp.iv.mpf) -> list[str]:
    return [mp.nstr(mp.mpf(value.a), 60), mp.nstr(mp.mpf(value.b), 60)]


def legacy_norm(x: mp.mpf, p: dict[str, mp.mpf]) -> mp.mpf:
    return mp.cos(p["legacy_omega"] * x + p["legacy_phi"]) / mp.cos(p["legacy_phi"]) / (1 + p["legacy_beta_tors"] * x)


def strict_norm(d: mp.mpf, p: dict[str, mp.mpf]) -> mp.mpf:
    return mp.cos(p["strict_omega"] * d + p["strict_phi"]) / mp.cos(p["strict_phi"]) / (1 + p["strict_beta"] * d ** p["strict_eta"])


def legacy_derivative_interval() -> mp.iv.mpf:
    mp.iv.dps = IV_DPS
    p = mp_params()
    omega = iv_mpf(p["legacy_omega"])
    phi = iv_mpf(p["legacy_phi"])
    beta = iv_mpf(p["legacy_beta_tors"])
    x = mp.iv.mpf([BRANCH_LEFT, BRANCH_RIGHT])
    theta = omega * x + phi
    numerator = -omega * mp.iv.sin(theta) * (1 + beta * x) - beta * mp.iv.cos(theta)
    denominator = mp.iv.cos(phi) * (1 + beta * x) ** 2
    return numerator / denominator


def branch_bracket_rows(d_values: list[str]) -> list[dict[str, Any]]:
    mp.mp.dps = POINT_DPS
    p = mp_params()
    left_value = legacy_norm(BRANCH_LEFT, p)
    right_value = legacy_norm(BRANCH_RIGHT, p)
    rows = []
    for d_text in d_values:
        strict_value = strict_norm(mp.mpf(d_text), p)
        left_margin = left_value - strict_value
        right_margin = strict_value - right_value
        rows.append({
            "d": d_text,
            "strict_norm_value": mp.nstr(strict_value, 60),
            "left_margin_L0_minus_S": mp.nstr(left_margin, 60),
            "right_margin_S_minus_L2": mp.nstr(right_margin, 60),
            "inside_legacy_branch_range": left_margin >= 0 and right_margin > 0,
            "endpoint_case_at_left": left_margin == 0,
        })
    return rows


def build_branch_certificate(p2495: dict[str, Any]) -> dict[str, Any]:
    derivative = legacy_derivative_interval()
    derivative_bounds = interval_bounds_text(derivative)
    sample_rows = branch_bracket_rows(SAMPLE_D_STRINGS)
    z12_rows = branch_bracket_rows(Z12_D_STRINGS)
    mp.mp.dps = POINT_DPS
    p = mp_params()
    left_value = legacy_norm(BRANCH_LEFT, p)
    right_value = legacy_norm(BRANCH_RIGHT, p)
    all_rows = sample_rows + z12_rows
    return {
        "branch_equation": "L_norm(x(d)) = S_strict_norm(d)",
        "legacy_branch_interval": [mp.nstr(BRANCH_LEFT, 20), mp.nstr(BRANCH_RIGHT, 20)],
        "legacy_norm_left_value": mp.nstr(left_value, 60),
        "legacy_norm_right_value": mp.nstr(right_value, 60),
        "legacy_derivative_interval_backend": "mpmath.iv",
        "legacy_derivative_interval_dps": IV_DPS,
        "legacy_derivative_interval_on_0_2": derivative_bounds,
        "legacy_derivative_interval_strictly_negative": mp.mpf(derivative.b) < 0,
        "sample_branch_bracket_rows": sample_rows,
        "z12_branch_bracket_rows": z12_rows,
        "all_sample_targets_inside_branch_range": all(row["inside_legacy_branch_range"] for row in sample_rows),
        "all_z12_targets_inside_branch_range": all(row["inside_legacy_branch_range"] for row in z12_rows),
        "positive_right_margin_min": mp.nstr(min(mp.mpf(row["right_margin_S_minus_L2"]) for row in all_rows), 60),
        "nonzero_left_margin_min_excluding_d0": mp.nstr(min(mp.mpf(row["left_margin_L0_minus_S"]) for row in all_rows if row["d"] != "0"), 60),
        "p2495_interval_enclosure_inherited": p2495.get("all_sample_curvature_intervals_strictly_signed") is True,
        "formal_directed_rounding_backend_exported": False,
        "global_branch_theorem_exported": False,
        "curvature_dynamic_source_exported": False,
        "legacy_to_strict_bridge_atom_exported": False,
        "strict_compression_dynamic_source_exported": False,
        "selector_source_theorem_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "toe_closure_exported": False,
        "certificate_fingerprint_sha256": sha256_json([derivative_bounds, sample_rows, z12_rows]),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2496/S1446 phase-normalized inverse-branch existence/uniqueness certificate

`P2496/S1446` audits the inverse branch used by P2493--P2495.  On the legacy normalized branch interval `x in [0,2]`, `mpmath.iv` encloses the legacy derivative in a strictly negative interval, so the audited branch is injective on that interval.  The strict normalized targets for the ten curvature samples and the twelve Z12 nodes lie in the corresponding legacy output range, with positive right-end margins and only the expected `d=0` left-end equality.  Thus the finite inverse-branch calls used by the curvature certificates have an existence/uniqueness bracket rather than being unbracketed point solves.

This is a finite interval/bracketing certificate, not a global analytic branch theorem or a formal directed-rounding proof backend.  It exports no curvature source, no damping bridge atom, no strict compression source, no selector/source closure, no QW-2191 discharge, no role-transfer, no role-bearing `L_total`, no physical-value generation, and no ToE closure.
"""
    lag_section = """
## P2496/S1446 inverse-branch existence/uniqueness guard

`P2496/S1446` adds a finite bracketing guard behind the P2493--P2495 curvature audits: the legacy normalized branch is interval-monotone on `[0,2]`, and all audited strict normalized targets are bracketed in its range.  This licenses the finite inverse-branch evaluations used by the curvature guards, but does not export a global branch theorem, nonlinear compression-flow source, bridge atom, role-transfer theorem, QW-2191 discharge, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2496/S1446 phase-normalized inverse-branch existence/uniqueness certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2496/S1446 inverse-branch existence/uniqueness guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2495 = theorem(sources["P2495_INTERVAL_ENCLOSURE"], "phase_normalized_compression_curvature_interval_enclosure_certificate")
    cert = build_branch_certificate(p2495)
    theorem_export = {
        "theorem_name": "P2496_T1_phase_normalized_inverse_branch_interval_existence_uniqueness_certificate",
        "audited_chain": ["P2493/S1443", "P2494/S1444", "P2495/S1445"],
        "inverse_branch_interval_existence_uniqueness_certificate": cert,
        "legacy_derivative_interval_strictly_negative": cert["legacy_derivative_interval_strictly_negative"],
        "all_sample_targets_inside_branch_range": cert["all_sample_targets_inside_branch_range"],
        "all_z12_targets_inside_branch_range": cert["all_z12_targets_inside_branch_range"],
        "positive_right_margin_min": cert["positive_right_margin_min"],
        "nonzero_left_margin_min_excluding_d0": cert["nonzero_left_margin_min_excluding_d0"],
        "p2495_interval_enclosure_inherited": cert["p2495_interval_enclosure_inherited"],
        "formal_directed_rounding_backend_exported": False,
        "global_branch_theorem_exported": False,
        "curvature_dynamic_source_exported": False,
        "legacy_to_strict_bridge_atom_exported": False,
        "strict_compression_dynamic_source_exported": False,
        "selector_source_theorem_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "P2496 brackets the finite inverse-branch evaluations used by P2493-P2495; it is not a global branch theorem.",
            "The interval derivative enclosure is numerical interval evidence, not a formal directed-rounding backend for real analysis.",
            "No damping bridge atom, strict compression dynamic source, selector/source theorem, QW-2191 discharge, role-transfer theorem, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Either promote the finite branch monotonicity/bracketing certificate to a formal directed-rounding theorem, or derive the nonlinear compression-flow ODE/source from strict-side geometry; branch bracketing alone does not close the bridge.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "legacy_branch_monotone_checked": theorem_export["legacy_derivative_interval_strictly_negative"],
        "sample_targets_bracketed": theorem_export["all_sample_targets_inside_branch_range"],
        "z12_targets_bracketed": theorem_export["all_z12_targets_inside_branch_range"],
        "margins_recorded": mp.mpf(theorem_export["positive_right_margin_min"]) > 0 and mp.mpf(theorem_export["nonzero_left_margin_min_excluding_d0"]) > 0,
        "p2495_interval_enclosure_inherited": theorem_export["p2495_interval_enclosure_inherited"],
        "no_closure_inflation": not any(theorem_export[key] for key in [
            "formal_directed_rounding_backend_exported",
            "global_branch_theorem_exported",
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
        "packet_id": "P2496",
        "stage_id": "S1446",
        "status": "PHASE_NORMALIZED_INVERSE_BRANCH_INTERVAL_EXISTENCE_UNIQUENESS_CERTIFICATE_NO_GLOBAL_BRANCH_THEOREM_NO_SOURCE_EXPORT_NO_BRIDGE_ATOM_NO_QW2191_NO_ROLE_TRANSFER_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "phase_normalized_inverse_branch_interval_existence_uniqueness_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["phase_normalized_inverse_branch_interval_existence_uniqueness_certificate"]["theorem_export"]
    cert = t["inverse_branch_interval_existence_uniqueness_certificate"]
    lines = [
        "# P2496/S1446 phase-normalized inverse-branch interval existence/uniqueness certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Branch certificate summary",
        "",
        f"Audited chain: `{', '.join(t['audited_chain'])}`.",
        f"Legacy branch interval: `{cert['legacy_branch_interval']}`.",
        f"Legacy derivative interval on `[0,2]`: `{cert['legacy_derivative_interval_on_0_2']}`.",
        f"Derivative interval strictly negative: `{t['legacy_derivative_interval_strictly_negative']}`.",
        f"All sample targets bracketed: `{t['all_sample_targets_inside_branch_range']}`.",
        f"All Z12 targets bracketed: `{t['all_z12_targets_inside_branch_range']}`.",
        f"Minimum positive right margin: `{t['positive_right_margin_min']}`.",
        f"Minimum nonzero left margin excluding `d=0`: `{t['nonzero_left_margin_min_excluding_d0']}`.",
        "",
        "## Sample branch rows",
        "",
    ]
    for row in cert["sample_branch_bracket_rows"]:
        lines.append(f"- `d={row['d']}`: inside `{row['inside_legacy_branch_range']}`, left margin `{row['left_margin_L0_minus_S']}`, right margin `{row['right_margin_S_minus_L2']}`.")
    lines += [
        "",
        "## Negative controls",
        "",
        "P2496 is not a global branch theorem or formal directed-rounding proof backend and does not export a curvature dynamic source, bridge atom, strict compression source theorem, selector/source theorem, QW-2191 discharge, role-transfer license, physical-value generator, or ToE closure.",
        "",
        "## Lay summary",
        "",
        "P2496 checks that the inverse branch used by P2493-P2495 is actually bracketed and unique on the audited finite domain.  It supports the curvature audit plumbing without deriving the strict compression source.",
        "",
        "## Fingerprints",
        "",
        f"Certificate fingerprint: `{cert['certificate_fingerprint_sha256']}`.",
        f"Theorem fingerprint: `{payload['phase_normalized_inverse_branch_interval_existence_uniqueness_certificate']['theorem_fingerprint_sha256']}`.",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
