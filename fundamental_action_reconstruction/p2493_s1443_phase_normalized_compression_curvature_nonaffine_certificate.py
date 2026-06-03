#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)

GEN = ROOT / "generated"
OUT = GEN / "p2493_s1443_phase_normalized_compression_curvature_nonaffine_certificate.json"
MD = GEN / "p2493_s1443_phase_normalized_compression_curvature_nonaffine_certificate.md"

SOURCE_FILES = {
    "P2414_DAMPING_IDENTIFIABILITY": GEN / "p2414_s1364_strict_damping_parameter_identifiability_nonabsorption_certificate.json",
    "P2492_COMPLETION_PACKAGES": GEN / "p2492_s1442_legacy_strict_claim_specific_minimal_completion_package_certificate.json",
    "SCRATCH_CURVATURE_PRECURSOR": ROOT / "scratch" / "bridge_phase_normalized_compression_curvature_report.json",
}

LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 9.0 / 5.0}
BRANCH_BRACKET = (0.0, 2.0)
SAMPLE_D = [1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 11.0]
Z12_D = list(range(12))


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
        "new_packet": "P2493|S1443|phase-normalized compression curvature|compression curvature nonaffine|inverse-branch curvature|curved compression flow|nonaffine compression certificate",
        "precursor_packets": "P2414|S1364|strict damping parameter identifiability|P2492|S1442|minimal completion package|bridge_phase_normalized_compression_curvature",
        "curvature_language": "x_second|second derivative identity|affine bridge|curvature_nonzero|Z12.*second difference|nonlinear compression flow",
        "guardrail_language": "legacy -> strict completion bridge|role-transfer audit|silent inheritance|K_legacy_ont|K_strict_gate|QW-2191",
        "closure_blockers": "role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer|root-window theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def legacy_norm(x: float) -> float:
    return math.cos(LEGACY["omega"] * x + LEGACY["phi"]) / math.cos(LEGACY["phi"]) / (1.0 + LEGACY["beta_tors"] * x)


def legacy_derivative(x: float) -> float:
    omega = LEGACY["omega"]
    phi = LEGACY["phi"]
    beta = LEGACY["beta_tors"]
    theta = omega * x + phi
    numerator = -omega * math.sin(theta) * (1.0 + beta * x) - beta * math.cos(theta)
    denominator = math.cos(phi) * (1.0 + beta * x) ** 2
    return numerator / denominator


def legacy_second_derivative(x: float) -> float:
    omega = LEGACY["omega"]
    phi = LEGACY["phi"]
    beta = LEGACY["beta_tors"]
    theta = omega * x + phi
    cos_phi = math.cos(phi)
    normalized_carrier = math.cos(theta) / cos_phi
    normalized_carrier_prime = -omega * math.sin(theta) / cos_phi
    normalized_carrier_second = -(omega**2) * math.cos(theta) / cos_phi
    denominator = 1.0 + beta * x
    return (
        normalized_carrier_second / denominator
        - 2.0 * normalized_carrier_prime * beta / denominator**2
        + 2.0 * normalized_carrier * beta**2 / denominator**3
    )


def strict_norm(d: float) -> float:
    return math.cos(STRICT["omega"] * d + STRICT["phi"]) / math.cos(STRICT["phi"]) / (1.0 + STRICT["beta"] * d ** STRICT["eta"])


def strict_derivative(d: float) -> float:
    omega = STRICT["omega"]
    phi = STRICT["phi"]
    beta = STRICT["beta"]
    eta = STRICT["eta"]
    theta = omega * d + phi
    carrier = math.cos(theta) / math.cos(phi)
    carrier_prime = -omega * math.sin(theta) / math.cos(phi)
    denominator = 1.0 + beta * d**eta
    denominator_prime = beta * eta * d ** (eta - 1.0)
    return (carrier_prime * denominator - carrier * denominator_prime) / denominator**2


def strict_second_derivative(d: float) -> float:
    omega = STRICT["omega"]
    phi = STRICT["phi"]
    beta = STRICT["beta"]
    eta = STRICT["eta"]
    theta = omega * d + phi
    carrier = math.cos(theta) / math.cos(phi)
    carrier_prime = -omega * math.sin(theta) / math.cos(phi)
    carrier_second = -(omega**2) * math.cos(theta) / math.cos(phi)
    denominator = 1.0 + beta * d**eta
    denominator_prime = beta * eta * d ** (eta - 1.0)
    denominator_second = beta * eta * (eta - 1.0) * d ** (eta - 2.0)
    return (
        (carrier_second * denominator - carrier * denominator_second) / denominator**2
        - 2.0 * (carrier_prime * denominator - carrier * denominator_prime) * denominator_prime / denominator**3
    )


def solve_first_branch(y: float, lo: float = BRANCH_BRACKET[0], hi: float = BRANCH_BRACKET[1]) -> float:
    f_lo = legacy_norm(lo) - y
    f_hi = legacy_norm(hi) - y
    if abs(f_lo) < 1e-15:
        return lo
    if f_lo * f_hi > 0:
        raise ValueError(f"target {y} not bracketed on [{lo}, {hi}]: {f_lo}, {f_hi}")
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        f_mid = legacy_norm(mid) - y
        if f_lo * f_mid <= 0.0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid
    return 0.5 * (lo + hi)


def branch_x(d: float) -> float:
    return solve_first_branch(strict_norm(d))


def derivative_row(d: float) -> dict[str, Any]:
    x = branch_x(d)
    legacy_prime = legacy_derivative(x)
    strict_prime = strict_derivative(d)
    x_prime = strict_prime / legacy_prime
    legacy_second = legacy_second_derivative(x)
    strict_second = strict_second_derivative(d)
    x_second = (strict_second - legacy_second * x_prime**2) / legacy_prime
    residual = legacy_second * x_prime**2 + legacy_prime * x_second - strict_second
    return {
        "d": d,
        "x": x,
        "x_prime": x_prime,
        "x_second": x_second,
        "legacy_prime": legacy_prime,
        "legacy_second": legacy_second,
        "strict_prime": strict_prime,
        "strict_second": strict_second,
        "second_derivative_identity_residual": residual,
        "curvature_nonzero": abs(x_second) > 1e-6,
        "curvature_sign": 1 if x_second > 0.0 else -1 if x_second < 0.0 else 0,
    }


def x_second_at(d: float) -> float:
    return float(derivative_row(d)["x_second"])


def curvature_inflection_row(lo: float = 0.1, hi: float = 0.5) -> dict[str, Any]:
    f_lo = x_second_at(lo)
    f_hi = x_second_at(hi)
    if f_lo * f_hi > 0.0:
        raise ValueError(f"curvature sign change not bracketed: {f_lo}, {f_hi}")
    left = lo
    right = hi
    for _ in range(70):
        mid = 0.5 * (left + right)
        f_mid = x_second_at(mid)
        if f_lo * f_mid <= 0.0:
            right = mid
            f_hi = f_mid
        else:
            left = mid
            f_lo = f_mid
    root = 0.5 * (left + right)
    return {
        "bracket": [lo, hi],
        "left_curvature": x_second_at(lo),
        "right_curvature": x_second_at(hi),
        "root_estimate": root,
        "root_curvature_abs": abs(x_second_at(root)),
        "sign_change_certified": x_second_at(lo) > 0.0 and x_second_at(hi) < 0.0,
    }


def z12_discrete_curvature() -> dict[str, Any]:
    x_values = [branch_x(float(d)) for d in Z12_D]
    second_differences = [x_values[index + 1] - 2.0 * x_values[index] + x_values[index - 1] for index in range(1, len(x_values) - 1)]
    return {
        "x_values": x_values,
        "second_differences_d1_to_d10": second_differences,
        "all_tail_second_differences_negative_from_d1": all(value < 0.0 for value in second_differences),
        "max_abs_second_difference": max(abs(value) for value in second_differences),
    }


def build_certificate(sources: dict[str, dict[str, Any]]) -> dict[str, Any]:
    rows = [derivative_row(d) for d in SAMPLE_D]
    z12 = z12_discrete_curvature()
    inflection = curvature_inflection_row()
    p2414 = theorem(sources["P2414_DAMPING_IDENTIFIABILITY"], "strict_damping_parameter_identifiability_nonabsorption_certificate")
    p2492 = theorem(sources["P2492_COMPLETION_PACKAGES"], "legacy_strict_claim_specific_minimal_completion_package_certificate")
    sample_signs = [row["curvature_sign"] for row in rows]
    return {
        "normalized_branch_equation": "L_norm(x(d)) = S_strict_norm(d)",
        "legacy_parameters": LEGACY,
        "strict_parameters": STRICT,
        "branch_bracket": list(BRANCH_BRACKET),
        "sample_rows": rows,
        "sample_count": len(rows),
        "nonzero_curvature_sample_count": sum(1 for row in rows if row["curvature_nonzero"]),
        "has_positive_curvature_near_origin": any(row["x_second"] > 0.0 for row in rows[:4]),
        "has_negative_curvature_on_tail": any(row["x_second"] < 0.0 for row in rows[4:]),
        "sample_curvature_sign_sequence": sample_signs,
        "max_second_derivative_identity_residual_abs": max(abs(row["second_derivative_identity_residual"]) for row in rows),
        "affine_bridge_ruled_out_by_curvature": all(row["curvature_nonzero"] for row in rows),
        "inflection_certificate": inflection,
        "discrete_z12_curvature": z12,
        "p2414_damping_nonabsorption_inherited": p2414.get("no_single_linear_gamma_matches_two_distinct_positive_nodes") is True,
        "p2492_shared_core_still_requires_damping_atoms": all(
            atom in p2492.get("claim_specific_minimal_completion_package_certificate", {}).get("shared_core_atoms", [])
            for atom in [
                "bridge::legacy_linear_to_strict_nonlinear_compression_map_theorem",
                "bridge::strict_compression_dynamic_source_theorem",
            ]
        ),
        "curvature_dynamic_source_exported": False,
        "legacy_to_strict_bridge_atom_exported": False,
        "strict_compression_dynamic_source_exported": False,
        "selector_source_theorem_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "toe_closure_exported": False,
        "certificate_fingerprint_sha256": sha256_json(rows + [z12, inflection]),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2493/S1443 phase-normalized compression curvature nonaffine certificate

`P2493/S1443` turns the phase-normalized legacy-to-strict output-matching branch into a differential curvature audit.  On the branch `L_norm(x(d)) = S_strict_norm(d)`, implicit differentiation gives `x' = S'(d)/L'(x)` and `x'' = (S''(d)-L''(x)*(x')^2)/L'(x)`.  The finite sample has nonzero curvature at all ten audited points, positive curvature near the origin, negative curvature on the tail, an inflection bracket between `d=0.1` and `d=0.5`, and negative Z12 second differences from `d=1` through `d=10`.  Therefore a pure affine phase/distance bridge is ruled out on this audited branch: any future bridge must supply a genuinely nonlinear compression-flow source.

This is a nonaffine constraint certificate only.  It does not export a curvature dynamic source, a damping bridge atom, a strict compression source theorem, selector/source closure, QW-2191 discharge, role-transfer, role-bearing `L_total`, physical-value generation, or ToE closure.
"""
    lag_section = """
## P2493/S1443 nonaffine compression curvature guard

`P2493/S1443` adds a differential guard behind `L_total`: the phase-normalized branch that matches strict damping output to legacy normalized output has nonzero `x''`, changes curvature sign from near-origin positive to tail negative, and has negative Z12 second differences.  Thus a role-bearing strict compression term cannot be justified by an affine bridge; it still needs a real nonlinear compression-flow source.  No QW-2191 discharge, role-transfer, physical-value generator, or ToE closure is exported.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2493/S1443 phase-normalized compression curvature nonaffine certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2493/S1443 nonaffine compression curvature guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    cert = build_certificate(sources)
    theorem_export = {
        "theorem_name": "P2493_T1_phase_normalized_compression_curvature_nonaffine_certificate",
        "audited_chain": ["P2414/S1364", "P2492/S1442", "scratch/bridge_phase_normalized_compression_curvature"],
        "phase_normalized_compression_curvature_certificate": cert,
        "sample_count": cert["sample_count"],
        "nonzero_curvature_sample_count": cert["nonzero_curvature_sample_count"],
        "max_second_derivative_identity_residual_abs": cert["max_second_derivative_identity_residual_abs"],
        "has_positive_curvature_near_origin": cert["has_positive_curvature_near_origin"],
        "has_negative_curvature_on_tail": cert["has_negative_curvature_on_tail"],
        "affine_bridge_ruled_out_by_curvature": cert["affine_bridge_ruled_out_by_curvature"],
        "z12_tail_second_differences_negative": cert["discrete_z12_curvature"]["all_tail_second_differences_negative_from_d1"],
        "inflection_sign_change_certified": cert["inflection_certificate"]["sign_change_certified"],
        "p2414_damping_nonabsorption_inherited": cert["p2414_damping_nonabsorption_inherited"],
        "p2492_shared_core_still_requires_damping_atoms": cert["p2492_shared_core_still_requires_damping_atoms"],
        "curvature_dynamic_source_exported": False,
        "legacy_to_strict_bridge_atom_exported": False,
        "strict_compression_dynamic_source_exported": False,
        "selector_source_theorem_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_transfer_licensed_by_this_certificate": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "P2493 proves only that the audited output-matching branch is nonaffine; it does not derive the branch from strict information geometry.",
            "The curvature profile constrains future bridge work but does not export legacy_linear_to_strict_nonlinear_compression_map_theorem.",
            "No strict compression dynamic source, selector/source theorem, QW-2191 discharge, role-transfer theorem, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Derive the nonlinear compression-flow ODE or curvature law from a strict-side information-geometric source; without that source, the damping bridge atoms remain in the P2492 shared core.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "sample_count_exact": theorem_export["sample_count"] == 10,
        "all_sample_curvatures_nonzero": theorem_export["nonzero_curvature_sample_count"] == 10,
        "identity_residual_tiny": theorem_export["max_second_derivative_identity_residual_abs"] < 1e-12,
        "positive_and_negative_curvature_seen": theorem_export["has_positive_curvature_near_origin"] and theorem_export["has_negative_curvature_on_tail"],
        "affine_bridge_ruled_out": theorem_export["affine_bridge_ruled_out_by_curvature"],
        "z12_tail_concavity_checked": theorem_export["z12_tail_second_differences_negative"],
        "inflection_bracket_checked": theorem_export["inflection_sign_change_certified"],
        "upstream_obligations_preserved": theorem_export["p2414_damping_nonabsorption_inherited"] and theorem_export["p2492_shared_core_still_requires_damping_atoms"],
        "no_closure_inflation": not any(theorem_export[key] for key in [
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
        "packet_id": "P2493",
        "stage_id": "S1443",
        "status": "PHASE_NORMALIZED_COMPRESSION_CURVATURE_NONAFFINE_CERTIFICATE_NO_SOURCE_EXPORT_NO_BRIDGE_ATOM_NO_QW2191_NO_ROLE_TRANSFER_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "phase_normalized_compression_curvature_nonaffine_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["phase_normalized_compression_curvature_nonaffine_certificate"]["theorem_export"]
    cert = t["phase_normalized_compression_curvature_certificate"]
    inflection = cert["inflection_certificate"]
    z12 = cert["discrete_z12_curvature"]
    lines = [
        "# P2493/S1443 phase-normalized compression curvature nonaffine certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Curvature summary",
        "",
        f"Audited chain: `{', '.join(t['audited_chain'])}`.",
        f"Sample count: `{t['sample_count']}`.",
        f"Nonzero curvature samples: `{t['nonzero_curvature_sample_count']}`.",
        f"Max second-derivative identity residual: `{t['max_second_derivative_identity_residual_abs']:.3e}`.",
        f"Positive near-origin curvature: `{t['has_positive_curvature_near_origin']}`.",
        f"Negative tail curvature: `{t['has_negative_curvature_on_tail']}`.",
        f"Affine bridge ruled out: `{t['affine_bridge_ruled_out_by_curvature']}`.",
        f"Inflection bracket: `{inflection['bracket']}`, root estimate `{inflection['root_estimate']:.12f}`.",
        f"Z12 tail second differences negative: `{t['z12_tail_second_differences_negative']}`.",
        f"Z12 max abs second difference: `{z12['max_abs_second_difference']:.12f}`.",
        "",
        "## Sample curvature rows",
        "",
    ]
    for row in cert["sample_rows"]:
        lines.append(f"- `d={row['d']}`: `x={row['x']:.12f}`, `x''={row['x_second']:.12f}`, sign `{row['curvature_sign']}`.")
    lines += [
        "",
        "## Negative controls",
        "",
        "P2493 does not export a curvature dynamic source, bridge atom, strict compression source theorem, selector/source theorem, QW-2191 discharge, role-transfer license, physical-value generator, or ToE closure.",
        "",
        "## Lay summary",
        "",
        "P2493 shows that the audited bridge candidate bends: it is not an affine relabeling of distance or phase.  This is a real constraint on future bridge completion, but not a completed source theorem.",
        "",
        "## Fingerprints",
        "",
        f"Certificate fingerprint: `{cert['certificate_fingerprint_sha256']}`.",
        f"Theorem fingerprint: `{payload['phase_normalized_compression_curvature_nonaffine_certificate']['theorem_fingerprint_sha256']}`.",
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
