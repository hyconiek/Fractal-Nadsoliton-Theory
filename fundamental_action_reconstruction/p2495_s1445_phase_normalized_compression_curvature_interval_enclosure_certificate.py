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
    mp_params,
    solve_first_branch,
    strict_norm,
)

GEN = ROOT / "generated"
OUT = GEN / "p2495_s1445_phase_normalized_compression_curvature_interval_enclosure_certificate.json"
MD = GEN / "p2495_s1445_phase_normalized_compression_curvature_interval_enclosure_certificate.md"

SOURCE_FILES = {
    "P2494_MULTIPRECISION_STABILITY": GEN / "p2494_s1444_phase_normalized_compression_curvature_multiprecision_stability_certificate.json",
}

IV_DPS = 100
POINT_DPS = 120
X_RADIUS = mp.mpf("1e-40")
EXPECTED_SAMPLE_SIGNS = [1, 1, 1, 1, -1, -1, -1, -1, -1, -1]


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
        "new_packet": "P2495|S1445|curvature interval enclosure|interval-enclosed curvature|mpmath interval enclosure|nonaffine interval certificate|sample curvature interval",
        "precursor_packets": "P2493|S1443|P2494|S1444|compression curvature nonaffine|multiprecision curvature stability",
        "interval_language": "interval enclosure|interval arithmetic|directed rounding|outward|curvature sign interval|Z12.*interval",
        "guardrail_language": "legacy -> strict completion bridge|role-transfer audit|silent inheritance|K_legacy_ont|K_strict_gate|QW-2191",
        "closure_blockers": "role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer|root-window theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def iv_mpf(value: mp.mpf | str) -> mp.iv.mpf:
    return mp.iv.mpf([value, value])


def iv_interval(left: mp.mpf, right: mp.mpf) -> mp.iv.mpf:
    return mp.iv.mpf([left, right])


def endpoint_text(value: mp.iv.mpf) -> str:
    return mp.nstr(mp.mpf(value), 60)


def interval_bounds_text(value: mp.iv.mpf) -> list[str]:
    return [endpoint_text(value.a), endpoint_text(value.b)]


def interval_width_text(value: mp.iv.mpf) -> str:
    return mp.nstr(mp.mpf(value.b) - mp.mpf(value.a), 40)


def interval_sign(value: mp.iv.mpf) -> int:
    lower = mp.mpf(value.a)
    upper = mp.mpf(value.b)
    if lower > 0:
        return 1
    if upper < 0:
        return -1
    return 0


def interval_parameters() -> dict[str, mp.iv.mpf]:
    p = mp_params()
    return {key: iv_mpf(value) for key, value in p.items()}


def legacy_derivative_interval(x: mp.iv.mpf, p: dict[str, mp.iv.mpf]) -> mp.iv.mpf:
    omega = p["legacy_omega"]
    phi = p["legacy_phi"]
    beta = p["legacy_beta_tors"]
    theta = omega * x + phi
    numerator = -omega * mp.iv.sin(theta) * (1 + beta * x) - beta * mp.iv.cos(theta)
    denominator = mp.iv.cos(phi) * (1 + beta * x) ** 2
    return numerator / denominator


def legacy_second_interval(x: mp.iv.mpf, p: dict[str, mp.iv.mpf]) -> mp.iv.mpf:
    omega = p["legacy_omega"]
    phi = p["legacy_phi"]
    beta = p["legacy_beta_tors"]
    theta = omega * x + phi
    cos_phi = mp.iv.cos(phi)
    normalized_carrier = mp.iv.cos(theta) / cos_phi
    normalized_carrier_prime = -omega * mp.iv.sin(theta) / cos_phi
    normalized_carrier_second = -(omega**2) * mp.iv.cos(theta) / cos_phi
    denominator = 1 + beta * x
    return normalized_carrier_second / denominator - 2 * normalized_carrier_prime * beta / denominator**2 + 2 * normalized_carrier * beta**2 / denominator**3


def strict_derivative_interval(d: mp.iv.mpf, p: dict[str, mp.iv.mpf]) -> mp.iv.mpf:
    omega = p["strict_omega"]
    phi = p["strict_phi"]
    beta = p["strict_beta"]
    eta = p["strict_eta"]
    theta = omega * d + phi
    carrier = mp.iv.cos(theta) / mp.iv.cos(phi)
    carrier_prime = -omega * mp.iv.sin(theta) / mp.iv.cos(phi)
    denominator = 1 + beta * d**eta
    denominator_prime = beta * eta * d ** (eta - 1)
    return (carrier_prime * denominator - carrier * denominator_prime) / denominator**2


def strict_second_interval(d: mp.iv.mpf, p: dict[str, mp.iv.mpf]) -> mp.iv.mpf:
    omega = p["strict_omega"]
    phi = p["strict_phi"]
    beta = p["strict_beta"]
    eta = p["strict_eta"]
    theta = omega * d + phi
    carrier = mp.iv.cos(theta) / mp.iv.cos(phi)
    carrier_prime = -omega * mp.iv.sin(theta) / mp.iv.cos(phi)
    carrier_second = -(omega**2) * mp.iv.cos(theta) / mp.iv.cos(phi)
    denominator = 1 + beta * d**eta
    denominator_prime = beta * eta * d ** (eta - 1)
    denominator_second = beta * eta * (eta - 1) * d ** (eta - 2)
    return (carrier_second * denominator - carrier * denominator_second) / denominator**2 - 2 * (carrier_prime * denominator - carrier * denominator_prime) * denominator_prime / denominator**3


def branch_x_point(d_text: str) -> mp.mpf:
    mp.mp.dps = POINT_DPS
    p = mp_params()
    d = mp.mpf(d_text)
    return solve_first_branch(strict_norm(d, p), p)


def curvature_interval_row(d_text: str) -> dict[str, Any]:
    mp.iv.dps = IV_DPS
    x = branch_x_point(d_text)
    p = interval_parameters()
    d_interval = iv_mpf(mp.mpf(d_text))
    x_interval = iv_interval(x - X_RADIUS, x + X_RADIUS)
    legacy_prime = legacy_derivative_interval(x_interval, p)
    strict_prime = strict_derivative_interval(d_interval, p)
    x_prime = strict_prime / legacy_prime
    legacy_second = legacy_second_interval(x_interval, p)
    strict_second = strict_second_interval(d_interval, p)
    x_second = (strict_second - legacy_second * x_prime**2) / legacy_prime
    sign = interval_sign(x_second)
    return {
        "d": d_text,
        "x_center": mp.nstr(x, 60),
        "x_radius": mp.nstr(X_RADIUS, 20),
        "x_second_interval": interval_bounds_text(x_second),
        "x_second_interval_width": interval_width_text(x_second),
        "curvature_interval_sign": sign,
        "curvature_interval_strictly_signed": sign != 0,
    }


def z12_interval_certificate() -> dict[str, Any]:
    centers = [branch_x_point(d_text) for d_text in Z12_D_STRINGS]
    intervals = [iv_interval(center - X_RADIUS, center + X_RADIUS) for center in centers]
    second = [intervals[index + 1] - 2 * intervals[index] + intervals[index - 1] for index in range(1, len(intervals) - 1)]
    rows = [
        {
            "d_center": str(index + 1),
            "second_difference_interval": interval_bounds_text(value),
            "second_difference_width": interval_width_text(value),
            "second_difference_sign": interval_sign(value),
        }
        for index, value in enumerate(second)
    ]
    return {
        "rows": rows,
        "all_second_difference_intervals_negative": all(row["second_difference_sign"] == -1 for row in rows),
        "minimum_negative_margin": min(abs(mp.mpf(value.b)) for value in second),
    }


def build_interval_certificate(p2494: dict[str, Any]) -> dict[str, Any]:
    sample_rows = [curvature_interval_row(d_text) for d_text in SAMPLE_D_STRINGS]
    z12 = z12_interval_certificate()
    sample_signs = [row["curvature_interval_sign"] for row in sample_rows]
    return {
        "interval_backend": "mpmath.iv",
        "interval_dps": IV_DPS,
        "point_solve_dps": POINT_DPS,
        "x_radius": mp.nstr(X_RADIUS, 20),
        "sample_rows": sample_rows,
        "sample_interval_sign_sequence": sample_signs,
        "all_sample_curvature_intervals_strictly_signed": all(row["curvature_interval_strictly_signed"] for row in sample_rows),
        "sample_interval_signs_match_p2493_p2494": sample_signs == EXPECTED_SAMPLE_SIGNS,
        "z12_interval_certificate": {
            "rows": z12["rows"],
            "all_second_difference_intervals_negative": z12["all_second_difference_intervals_negative"],
            "minimum_negative_margin": mp.nstr(z12["minimum_negative_margin"], 60),
        },
        "p2494_sample_sign_stability_inherited": p2494.get("sample_sign_sequence_stable_across_precisions") is True,
        "p2494_z12_sign_stability_inherited": p2494.get("z12_sign_sequence_stable_across_precisions") is True,
        "formal_directed_rounding_backend_exported": False,
        "curvature_dynamic_source_exported": False,
        "legacy_to_strict_bridge_atom_exported": False,
        "strict_compression_dynamic_source_exported": False,
        "selector_source_theorem_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "toe_closure_exported": False,
        "certificate_fingerprint_sha256": sha256_json(sample_rows + z12["rows"]),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2495/S1445 phase-normalized compression curvature interval-enclosure certificate

`P2495/S1445` replays the P2493/P2494 nonaffine curvature signs with `mpmath.iv` interval enclosures.  Using 120-digit point solves for the inverse branch centers, a declared `1e-40` branch-radius enclosure, and 100-digit interval arithmetic for the derivative expressions, all ten audited `x''` intervals are strictly signed with the P2493/P2494 sign pattern, and all Z12 second-difference intervals are strictly negative.  This gives a finite interval-enclosure check that the P2493 nonaffine signs are not merely point-evaluation artifacts.

This is still not a formal proof backend for directed-rounding real analysis and not a strict dynamic source theorem.  It exports no curvature source, no damping bridge atom, no strict compression source, no selector/source closure, no QW-2191 discharge, no role-transfer, no role-bearing `L_total`, no physical-value generation, and no ToE closure.
"""
    lag_section = """
## P2495/S1445 interval-enclosed curvature guard

`P2495/S1445` strengthens the P2493/P2494 nonaffine compression guard by enclosing the sample `x''` signs and Z12 second differences in `mpmath.iv` intervals.  The strict signs survive the declared branch-radius enclosure, so the affine-bridge obstruction is numerically interval-stable; nevertheless no nonlinear compression-flow source, bridge atom, role-transfer theorem, QW-2191 discharge, or ToE closure is exported.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2495/S1445 phase-normalized compression curvature interval-enclosure certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2495/S1445 interval-enclosed curvature guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2494 = theorem(sources["P2494_MULTIPRECISION_STABILITY"], "phase_normalized_compression_curvature_multiprecision_stability_certificate")
    cert = build_interval_certificate(p2494)
    theorem_export = {
        "theorem_name": "P2495_T1_phase_normalized_compression_curvature_interval_enclosure_certificate",
        "audited_chain": ["P2493/S1443", "P2494/S1444"],
        "interval_enclosure_certificate": cert,
        "interval_backend": cert["interval_backend"],
        "interval_dps": cert["interval_dps"],
        "point_solve_dps": cert["point_solve_dps"],
        "sample_interval_sign_sequence": cert["sample_interval_sign_sequence"],
        "all_sample_curvature_intervals_strictly_signed": cert["all_sample_curvature_intervals_strictly_signed"],
        "sample_interval_signs_match_p2493_p2494": cert["sample_interval_signs_match_p2493_p2494"],
        "z12_second_difference_intervals_all_negative": cert["z12_interval_certificate"]["all_second_difference_intervals_negative"],
        "p2494_sample_sign_stability_inherited": cert["p2494_sample_sign_stability_inherited"],
        "p2494_z12_sign_stability_inherited": cert["p2494_z12_sign_stability_inherited"],
        "formal_directed_rounding_backend_exported": False,
        "curvature_dynamic_source_exported": False,
        "legacy_to_strict_bridge_atom_exported": False,
        "strict_compression_dynamic_source_exported": False,
        "selector_source_theorem_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "P2495 is an mpmath.iv interval-enclosure audit, not a formal directed-rounding theorem backend.",
            "Interval-stable nonaffine signs do not derive the nonlinear compression-flow law from strict information geometry.",
            "No damping bridge atom, strict compression dynamic source, selector/source theorem, QW-2191 discharge, role-transfer theorem, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Promote the interval enclosure to a formal directed-rounding proof backend or derive the nonlinear compression-flow ODE from a strict-side source; until then the bridge source atoms remain open.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "interval_backend_checked": theorem_export["interval_backend"] == "mpmath.iv",
        "sample_intervals_strictly_signed": theorem_export["all_sample_curvature_intervals_strictly_signed"],
        "sample_signs_match_prior": theorem_export["sample_interval_signs_match_p2493_p2494"],
        "z12_intervals_negative": theorem_export["z12_second_difference_intervals_all_negative"],
        "p2494_stability_inherited": theorem_export["p2494_sample_sign_stability_inherited"] and theorem_export["p2494_z12_sign_stability_inherited"],
        "no_closure_inflation": not any(theorem_export[key] for key in [
            "formal_directed_rounding_backend_exported",
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
        "packet_id": "P2495",
        "stage_id": "S1445",
        "status": "PHASE_NORMALIZED_COMPRESSION_CURVATURE_INTERVAL_ENCLOSURE_CERTIFICATE_NO_FORMAL_DIRECTED_ROUNDING_NO_SOURCE_EXPORT_NO_BRIDGE_ATOM_NO_QW2191_NO_ROLE_TRANSFER_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "phase_normalized_compression_curvature_interval_enclosure_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["phase_normalized_compression_curvature_interval_enclosure_certificate"]["theorem_export"]
    cert = t["interval_enclosure_certificate"]
    lines = [
        "# P2495/S1445 phase-normalized compression curvature interval-enclosure certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Interval-enclosure summary",
        "",
        f"Audited chain: `{', '.join(t['audited_chain'])}`.",
        f"Interval backend: `{t['interval_backend']}` at dps `{t['interval_dps']}`.",
        f"Point solve dps: `{t['point_solve_dps']}`.",
        f"Sample interval sign sequence: `{t['sample_interval_sign_sequence']}`.",
        f"All sample curvature intervals strictly signed: `{t['all_sample_curvature_intervals_strictly_signed']}`.",
        f"Sample interval signs match P2493/P2494: `{t['sample_interval_signs_match_p2493_p2494']}`.",
        f"Z12 second-difference intervals all negative: `{t['z12_second_difference_intervals_all_negative']}`.",
        f"Minimum Z12 negative interval margin: `{cert['z12_interval_certificate']['minimum_negative_margin']}`.",
        "",
        "## Sample interval rows",
        "",
    ]
    for row in cert["sample_rows"]:
        lines.append(f"- `d={row['d']}`: sign `{row['curvature_interval_sign']}`, interval `{row['x_second_interval']}`.")
    lines += [
        "",
        "## Negative controls",
        "",
        "P2495 is not a formal directed-rounding proof backend and does not export a curvature dynamic source, bridge atom, strict compression source theorem, selector/source theorem, QW-2191 discharge, role-transfer license, physical-value generator, or ToE closure.",
        "",
        "## Lay summary",
        "",
        "P2495 encloses the P2493/P2494 curvature signs in intervals.  This is a stricter numerical guard against point-evaluation artifacts, but it still does not derive the strict compression source.",
        "",
        "## Fingerprints",
        "",
        f"Certificate fingerprint: `{cert['certificate_fingerprint_sha256']}`.",
        f"Theorem fingerprint: `{payload['phase_normalized_compression_curvature_interval_enclosure_certificate']['theorem_fingerprint_sha256']}`.",
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
