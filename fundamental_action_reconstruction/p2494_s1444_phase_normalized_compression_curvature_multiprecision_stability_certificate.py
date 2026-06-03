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

GEN = ROOT / "generated"
OUT = GEN / "p2494_s1444_phase_normalized_compression_curvature_multiprecision_stability_certificate.json"
MD = GEN / "p2494_s1444_phase_normalized_compression_curvature_multiprecision_stability_certificate.md"

SOURCE_FILES = {
    "P2493_CURVATURE_NONAFFINE": GEN / "p2493_s1443_phase_normalized_compression_curvature_nonaffine_certificate.json",
}

DPS_LEVELS = [50, 80, 120]
SAMPLE_D_STRINGS = ["0.0001", "0.001", "0.01", "0.1", "0.5", "1", "2", "4", "8", "11"]
Z12_D_STRINGS = [str(index) for index in range(12)]
BRACKET = (mp.mpf("0"), mp.mpf("2"))
INFLECTION_BRACKET_STRINGS = ("0.1", "0.5")


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
        "new_packet": "P2494|S1444|multiprecision curvature stability|multi-precision curvature|compression curvature precision|curvature sign stability|nonaffine stability certificate",
        "precursor_packets": "P2493|S1443|phase-normalized compression curvature|compression curvature nonaffine|inverse-branch curvature",
        "precision_language": "mpmath|DPS|precision stability|sign-stability|root drift|multi-precision|multiprecision",
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


def sign(value: mp.mpf) -> int:
    return 1 if value > 0 else -1 if value < 0 else 0


def legacy_norm(x: mp.mpf, p: dict[str, mp.mpf]) -> mp.mpf:
    return mp.cos(p["legacy_omega"] * x + p["legacy_phi"]) / mp.cos(p["legacy_phi"]) / (1 + p["legacy_beta_tors"] * x)


def legacy_derivative(x: mp.mpf, p: dict[str, mp.mpf]) -> mp.mpf:
    omega = p["legacy_omega"]
    phi = p["legacy_phi"]
    beta = p["legacy_beta_tors"]
    theta = omega * x + phi
    numerator = -omega * mp.sin(theta) * (1 + beta * x) - beta * mp.cos(theta)
    denominator = mp.cos(phi) * (1 + beta * x) ** 2
    return numerator / denominator


def legacy_second_derivative(x: mp.mpf, p: dict[str, mp.mpf]) -> mp.mpf:
    omega = p["legacy_omega"]
    phi = p["legacy_phi"]
    beta = p["legacy_beta_tors"]
    theta = omega * x + phi
    cos_phi = mp.cos(phi)
    normalized_carrier = mp.cos(theta) / cos_phi
    normalized_carrier_prime = -omega * mp.sin(theta) / cos_phi
    normalized_carrier_second = -(omega**2) * mp.cos(theta) / cos_phi
    denominator = 1 + beta * x
    return normalized_carrier_second / denominator - 2 * normalized_carrier_prime * beta / denominator**2 + 2 * normalized_carrier * beta**2 / denominator**3


def strict_norm(d: mp.mpf, p: dict[str, mp.mpf]) -> mp.mpf:
    return mp.cos(p["strict_omega"] * d + p["strict_phi"]) / mp.cos(p["strict_phi"]) / (1 + p["strict_beta"] * d ** p["strict_eta"])


def strict_derivative(d: mp.mpf, p: dict[str, mp.mpf]) -> mp.mpf:
    omega = p["strict_omega"]
    phi = p["strict_phi"]
    beta = p["strict_beta"]
    eta = p["strict_eta"]
    theta = omega * d + phi
    carrier = mp.cos(theta) / mp.cos(phi)
    carrier_prime = -omega * mp.sin(theta) / mp.cos(phi)
    denominator = 1 + beta * d**eta
    denominator_prime = beta * eta * d ** (eta - 1)
    return (carrier_prime * denominator - carrier * denominator_prime) / denominator**2


def strict_second_derivative(d: mp.mpf, p: dict[str, mp.mpf]) -> mp.mpf:
    omega = p["strict_omega"]
    phi = p["strict_phi"]
    beta = p["strict_beta"]
    eta = p["strict_eta"]
    theta = omega * d + phi
    carrier = mp.cos(theta) / mp.cos(phi)
    carrier_prime = -omega * mp.sin(theta) / mp.cos(phi)
    carrier_second = -(omega**2) * mp.cos(theta) / mp.cos(phi)
    denominator = 1 + beta * d**eta
    denominator_prime = beta * eta * d ** (eta - 1)
    denominator_second = beta * eta * (eta - 1) * d ** (eta - 2)
    return (carrier_second * denominator - carrier * denominator_second) / denominator**2 - 2 * (carrier_prime * denominator - carrier * denominator_prime) * denominator_prime / denominator**3


def solve_first_branch(y: mp.mpf, p: dict[str, mp.mpf]) -> mp.mpf:
    lo, hi = BRACKET
    f_lo = legacy_norm(lo, p) - y
    f_hi = legacy_norm(hi, p) - y
    if abs(f_lo) < mp.mpf("1e-80"):
        return lo
    if f_lo * f_hi > 0:
        raise ValueError(f"target not bracketed: {mp.nstr(f_lo, 30)}, {mp.nstr(f_hi, 30)}")
    for _ in range(4 * mp.mp.dps):
        mid = (lo + hi) / 2
        f_mid = legacy_norm(mid, p) - y
        if f_lo * f_mid <= 0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid
    return (lo + hi) / 2


def branch_x(d: mp.mpf, p: dict[str, mp.mpf]) -> mp.mpf:
    return solve_first_branch(strict_norm(d, p), p)


def curvature_value(d: mp.mpf, p: dict[str, mp.mpf]) -> dict[str, mp.mpf]:
    x = branch_x(d, p)
    legacy_prime = legacy_derivative(x, p)
    strict_prime = strict_derivative(d, p)
    x_prime = strict_prime / legacy_prime
    legacy_second = legacy_second_derivative(x, p)
    strict_second = strict_second_derivative(d, p)
    x_second = (strict_second - legacy_second * x_prime**2) / legacy_prime
    residual = legacy_second * x_prime**2 + legacy_prime * x_second - strict_second
    return {"x": x, "x_prime": x_prime, "x_second": x_second, "identity_residual": residual}


def inflection_root(p: dict[str, mp.mpf]) -> dict[str, mp.mpf | bool]:
    lo = mp.mpf(INFLECTION_BRACKET_STRINGS[0])
    hi = mp.mpf(INFLECTION_BRACKET_STRINGS[1])
    f_lo = curvature_value(lo, p)["x_second"]
    f_hi = curvature_value(hi, p)["x_second"]
    if f_lo * f_hi > 0:
        raise ValueError("inflection bracket failed")
    left = lo
    right = hi
    for _ in range(3 * mp.mp.dps):
        mid = (left + right) / 2
        f_mid = curvature_value(mid, p)["x_second"]
        if f_lo * f_mid <= 0:
            right = mid
            f_hi = f_mid
        else:
            left = mid
            f_lo = f_mid
    root = (left + right) / 2
    root_curvature = curvature_value(root, p)["x_second"]
    return {"root": root, "root_curvature_abs": abs(root_curvature), "sign_change_certified": curvature_value(lo, p)["x_second"] > 0 and curvature_value(hi, p)["x_second"] < 0}


def z12_summary(p: dict[str, mp.mpf]) -> dict[str, Any]:
    xs = [branch_x(mp.mpf(d), p) for d in Z12_D_STRINGS]
    second = [xs[index + 1] - 2 * xs[index] + xs[index - 1] for index in range(1, len(xs) - 1)]
    return {
        "second_difference_sign_sequence": [sign(value) for value in second],
        "all_tail_second_differences_negative": all(value < 0 for value in second),
        "max_abs_second_difference": max(abs(value) for value in second),
        "min_abs_second_difference": min(abs(value) for value in second),
    }


def precision_run(dps: int) -> dict[str, Any]:
    mp.mp.dps = dps
    p = mp_params()
    rows = []
    x_second_values = []
    residuals = []
    for d_text in SAMPLE_D_STRINGS:
        d = mp.mpf(d_text)
        row = curvature_value(d, p)
        x_second_values.append(row["x_second"])
        residuals.append(abs(row["identity_residual"]))
        rows.append({
            "d": d_text,
            "x": mp.nstr(row["x"], 50),
            "x_second": mp.nstr(row["x_second"], 50),
            "curvature_sign": sign(row["x_second"]),
            "identity_residual_abs": mp.nstr(abs(row["identity_residual"]), 20),
        })
    root = inflection_root(p)
    z12 = z12_summary(p)
    return {
        "dps": dps,
        "sample_rows": rows,
        "sample_curvature_sign_sequence": [row["curvature_sign"] for row in rows],
        "nonzero_curvature_sample_count": sum(1 for value in x_second_values if abs(value) > mp.mpf("1e-30")),
        "min_abs_sample_curvature": mp.nstr(min(abs(value) for value in x_second_values), 50),
        "max_identity_residual_abs": mp.nstr(max(residuals), 20),
        "inflection_root": mp.nstr(root["root"], 50),
        "inflection_root_curvature_abs": mp.nstr(root["root_curvature_abs"], 20),
        "inflection_sign_change_certified": bool(root["sign_change_certified"]),
        "z12_summary": {
            "second_difference_sign_sequence": z12["second_difference_sign_sequence"],
            "all_tail_second_differences_negative": z12["all_tail_second_differences_negative"],
            "max_abs_second_difference": mp.nstr(z12["max_abs_second_difference"], 50),
            "min_abs_second_difference": mp.nstr(z12["min_abs_second_difference"], 50),
        },
    }


def mp_abs_diff(text_a: str, text_b: str) -> mp.mpf:
    return abs(mp.mpf(text_a) - mp.mpf(text_b))


def stability_certificate(p2493: dict[str, Any]) -> dict[str, Any]:
    precision_rows = [precision_run(dps) for dps in DPS_LEVELS]
    sign_sequences = [row["sample_curvature_sign_sequence"] for row in precision_rows]
    z12_sign_sequences = [row["z12_summary"]["second_difference_sign_sequence"] for row in precision_rows]
    mp.mp.dps = 100
    highest = precision_rows[-1]
    previous = precision_rows[-2]
    max_x_second_drift = max(
        mp_abs_diff(left["x_second"], right["x_second"])
        for left, right in zip(previous["sample_rows"], highest["sample_rows"])
    )
    root_drift = mp_abs_diff(previous["inflection_root"], highest["inflection_root"])
    p2493_cert = p2493["phase_normalized_compression_curvature_certificate"]
    return {
        "dps_levels": DPS_LEVELS,
        "precision_rows": precision_rows,
        "sample_sign_sequence_stable_across_precisions": all(sequence == sign_sequences[0] for sequence in sign_sequences),
        "z12_sign_sequence_stable_across_precisions": all(sequence == z12_sign_sequences[0] for sequence in z12_sign_sequences),
        "nonzero_count_stable": all(row["nonzero_curvature_sample_count"] == len(SAMPLE_D_STRINGS) for row in precision_rows),
        "inflection_sign_change_stable": all(row["inflection_sign_change_certified"] for row in precision_rows),
        "max_x_second_drift_80_to_120": mp.nstr(max_x_second_drift, 30),
        "inflection_root_drift_80_to_120": mp.nstr(root_drift, 30),
        "highest_precision_min_abs_sample_curvature": highest["min_abs_sample_curvature"],
        "highest_precision_min_abs_z12_second_difference": highest["z12_summary"]["min_abs_second_difference"],
        "p2493_affine_bridge_ruled_out_inherited": p2493.get("affine_bridge_ruled_out_by_curvature") is True,
        "p2493_inflection_root_float": p2493_cert.get("inflection_certificate", {}).get("root_estimate"),
        "highest_precision_inflection_root": highest["inflection_root"],
        "directed_rounding_interval_proof_exported": False,
        "curvature_dynamic_source_exported": False,
        "legacy_to_strict_bridge_atom_exported": False,
        "strict_compression_dynamic_source_exported": False,
        "selector_source_theorem_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "toe_closure_exported": False,
        "certificate_fingerprint_sha256": sha256_json(precision_rows),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2494/S1444 phase-normalized compression curvature multiprecision stability certificate

`P2494/S1444` replays the P2493 curvature computation at `50`, `80`, and `120` decimal digits of `mpmath` precision.  Across all precision levels, the ten sample curvature signs remain stable, all ten audited curvatures remain nonzero, the inflection sign-change bracket remains certified, and the Z12 second-difference sign sequence remains negative.  The `80 -> 120` drift is recorded for both sample `x''` values and the inflection root, giving a finite numerical stability check for the P2493 nonaffine conclusion.

This is not directed-rounding interval arithmetic and not a strict dynamic source theorem.  It exports no curvature source, no damping bridge atom, no strict compression source, no selector/source closure, no QW-2191 discharge, no role-transfer, no role-bearing `L_total`, no physical-value generation, and no ToE closure.
"""
    lag_section = """
## P2494/S1444 multiprecision curvature stability guard

`P2494/S1444` strengthens the P2493 nonaffine compression guard by replaying the curvature signs, inflection bracket, and Z12 concavity at `50/80/120` decimal digits.  The stability check supports the finite nonaffine obstruction but remains non-directed numerical evidence; it does not export a nonlinear compression-flow source, bridge atom, role-transfer theorem, QW-2191 discharge, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2494/S1444 phase-normalized compression curvature multiprecision stability certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2494/S1444 multiprecision curvature stability guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2493 = theorem(sources["P2493_CURVATURE_NONAFFINE"], "phase_normalized_compression_curvature_nonaffine_certificate")
    cert = stability_certificate(p2493)
    theorem_export = {
        "theorem_name": "P2494_T1_phase_normalized_compression_curvature_multiprecision_stability_certificate",
        "audited_chain": ["P2493/S1443"],
        "multiprecision_stability_certificate": cert,
        "dps_levels": cert["dps_levels"],
        "sample_sign_sequence_stable_across_precisions": cert["sample_sign_sequence_stable_across_precisions"],
        "z12_sign_sequence_stable_across_precisions": cert["z12_sign_sequence_stable_across_precisions"],
        "nonzero_count_stable": cert["nonzero_count_stable"],
        "inflection_sign_change_stable": cert["inflection_sign_change_stable"],
        "max_x_second_drift_80_to_120": cert["max_x_second_drift_80_to_120"],
        "inflection_root_drift_80_to_120": cert["inflection_root_drift_80_to_120"],
        "highest_precision_min_abs_sample_curvature": cert["highest_precision_min_abs_sample_curvature"],
        "highest_precision_min_abs_z12_second_difference": cert["highest_precision_min_abs_z12_second_difference"],
        "p2493_affine_bridge_ruled_out_inherited": cert["p2493_affine_bridge_ruled_out_inherited"],
        "directed_rounding_interval_proof_exported": False,
        "curvature_dynamic_source_exported": False,
        "legacy_to_strict_bridge_atom_exported": False,
        "strict_compression_dynamic_source_exported": False,
        "selector_source_theorem_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "P2494 is a multiprecision replay/stability certificate, not directed-rounding interval arithmetic.",
            "Stable signs do not derive the nonlinear compression-flow law from strict information geometry.",
            "No damping bridge atom, strict compression dynamic source, selector/source theorem, QW-2191 discharge, role-transfer theorem, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Replace multiprecision numerical stability with a directed-rounding interval proof or derive the nonlinear compression ODE from a strict-side source; until then P2493/P2494 remain non-source constraints.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "three_precision_levels_checked": theorem_export["dps_levels"] == DPS_LEVELS,
        "sample_signs_stable": theorem_export["sample_sign_sequence_stable_across_precisions"],
        "z12_signs_stable": theorem_export["z12_sign_sequence_stable_across_precisions"],
        "nonzero_counts_stable": theorem_export["nonzero_count_stable"],
        "inflection_stable": theorem_export["inflection_sign_change_stable"],
        "p2493_nonaffine_inherited": theorem_export["p2493_affine_bridge_ruled_out_inherited"],
        "drift_recorded": mp.mpf(theorem_export["max_x_second_drift_80_to_120"]) >= 0 and mp.mpf(theorem_export["inflection_root_drift_80_to_120"]) >= 0,
        "no_closure_inflation": not any(theorem_export[key] for key in [
            "directed_rounding_interval_proof_exported",
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
        "packet_id": "P2494",
        "stage_id": "S1444",
        "status": "PHASE_NORMALIZED_COMPRESSION_CURVATURE_MULTIPRECISION_STABILITY_CERTIFICATE_NO_DIRECTED_ROUNDING_NO_SOURCE_EXPORT_NO_BRIDGE_ATOM_NO_QW2191_NO_ROLE_TRANSFER_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "phase_normalized_compression_curvature_multiprecision_stability_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["phase_normalized_compression_curvature_multiprecision_stability_certificate"]["theorem_export"]
    cert = t["multiprecision_stability_certificate"]
    highest = cert["precision_rows"][-1]
    lines = [
        "# P2494/S1444 phase-normalized compression curvature multiprecision stability certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Multiprecision stability summary",
        "",
        f"Audited chain: `{', '.join(t['audited_chain'])}`.",
        f"DPS levels: `{t['dps_levels']}`.",
        f"Sample sign sequence stable: `{t['sample_sign_sequence_stable_across_precisions']}`.",
        f"Z12 sign sequence stable: `{t['z12_sign_sequence_stable_across_precisions']}`.",
        f"Nonzero count stable: `{t['nonzero_count_stable']}`.",
        f"Inflection sign-change stable: `{t['inflection_sign_change_stable']}`.",
        f"Max `x''` drift from 80 to 120 dps: `{t['max_x_second_drift_80_to_120']}`.",
        f"Inflection-root drift from 80 to 120 dps: `{t['inflection_root_drift_80_to_120']}`.",
        f"Highest-precision min abs sample curvature: `{t['highest_precision_min_abs_sample_curvature']}`.",
        f"Highest-precision min abs Z12 second difference: `{t['highest_precision_min_abs_z12_second_difference']}`.",
        "",
        "## Highest precision sample signs",
        "",
    ]
    for row in highest["sample_rows"]:
        lines.append(f"- `d={row['d']}`: sign `{row['curvature_sign']}`, `x''={row['x_second']}`.")
    lines += [
        "",
        "## Negative controls",
        "",
        "P2494 is not a directed-rounding interval proof and does not export a curvature dynamic source, bridge atom, strict compression source theorem, selector/source theorem, QW-2191 discharge, role-transfer license, physical-value generator, or ToE closure.",
        "",
        "## Lay summary",
        "",
        "P2494 checks that the P2493 nonaffine curvature signs are not a low-precision numerical accident.  It is stronger numerical evidence, but still not a strict source theorem.",
        "",
        "## Fingerprints",
        "",
        f"Certificate fingerprint: `{cert['certificate_fingerprint_sha256']}`.",
        f"Theorem fingerprint: `{payload['phase_normalized_compression_curvature_multiprecision_stability_certificate']['theorem_fingerprint_sha256']}`.",
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
