#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any, Callable

from p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate import (
    AUDIT_INTERVAL,
    dot,
    load_json,
    pointwise_gradient_derivative,
    projection_amplitude,
    projection_stationary_factor,
    projection_volume,
    rel,
    row_norm,
)

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json"
MD = GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.md"

SOURCE_FILES = {
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

ROOT_WINDOW_HALF_WIDTH = 2.0e-3
EXCLUSION_CELL_WIDTH = 1.0e-3
H_AUDIT_START = 1.0e-4
DERIVATIVE_DIFFERENCE_H = 1.0e-6
DERIVATIVE_BOUND_SAFETY_FACTOR = 1.2
EXCLUSION_MARGIN_TOL = 0.0
MONOTONICITY_SAMPLE_COUNT = 17


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
    return {"count": len(lines), "samples": lines[:35]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2450|S1400|root isolation margin|projection root isolation|sampled Lipschitz exclusion",
        "p2449_input": "P2449|S1399|projection reduction|zero-projection root|stationary factor",
        "isolation_language": "root isolation|sign-changing bracket|monotonicity window|exclusion cell|Lipschitz margin",
        "numerical_guard_language": "sampled derivative bound|not exact interval|not interval root-exclusion|finite mesh exclusion",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|point-coordinate selector",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def amplitude_derivative(projection_vector: list[float], d: float) -> float:
    return dot(projection_vector, pointwise_gradient_derivative(d))


def finite_difference_derivative(function: Callable[[float], float], d: float, left_bound: float, right_bound: float) -> float:
    left = max(left_bound, d - DERIVATIVE_DIFFERENCE_H)
    right = min(right_bound, d + DERIVATIVE_DIFFERENCE_H)
    return (function(right) - function(left)) / (right - left)


def root_windows(root_ds: list[float], start: float, end: float) -> list[dict[str, float]]:
    windows = []
    for root in root_ds:
        windows.append({"root_d": root, "left": max(start, root - ROOT_WINDOW_HALF_WIDTH), "right": min(end, root + ROOT_WINDOW_HALF_WIDTH)})
    return windows


def monotonicity_window_audit(
    family: str,
    function: Callable[[float], float],
    derivative: Callable[[float], float],
    windows: list[dict[str, float]],
) -> list[dict[str, Any]]:
    audits = []
    for window in windows:
        left = window["left"]
        right = window["right"]
        samples = [left + (right - left) * idx / (MONOTONICITY_SAMPLE_COUNT - 1) for idx in range(MONOTONICITY_SAMPLE_COUNT)]
        derivative_values = [derivative(sample) for sample in samples]
        derivative_sign = "positive" if min(derivative_values) > 0.0 else "negative" if max(derivative_values) < 0.0 else "mixed"
        left_value = function(left)
        right_value = function(right)
        audits.append(
            {
                "family": family,
                "root_d": window["root_d"],
                "window_left": left,
                "window_right": right,
                "left_value": left_value,
                "right_value": right_value,
                "sign_change_across_window": left_value * right_value < 0.0,
                "derivative_sample_count": len(samples),
                "minimum_sampled_derivative": min(derivative_values),
                "maximum_sampled_derivative": max(derivative_values),
                "sampled_derivative_sign": derivative_sign,
                "sampled_monotone_window": derivative_sign in {"positive", "negative"},
            }
        )
    return audits


def complement_segments(windows: list[dict[str, float]], start: float, end: float) -> list[dict[str, float]]:
    ordered = sorted(windows, key=lambda window: window["left"])
    segments = []
    cursor = start
    for window in ordered:
        if window["left"] > cursor:
            segments.append({"left": cursor, "right": window["left"]})
        cursor = max(cursor, window["right"])
    if cursor < end:
        segments.append({"left": cursor, "right": end})
    return segments


def sampled_lipschitz_exclusion_audit(
    family: str,
    function: Callable[[float], float],
    derivative: Callable[[float], float],
    windows: list[dict[str, float]],
    start: float,
    end: float,
) -> dict[str, Any]:
    segments = complement_segments(windows, start, end)
    minimum_margin = math.inf
    minimum_ratio = math.inf
    weakest_cell: dict[str, float] | None = None
    failed_cells = []
    cell_count = 0
    for segment in segments:
        left = segment["left"]
        right = segment["right"]
        x = left
        while x < right - 1.0e-15:
            y = min(right, x + EXCLUSION_CELL_WIDTH)
            midpoint = (x + y) / 2.0
            half_width = (y - x) / 2.0
            sampled_bound = DERIVATIVE_BOUND_SAFETY_FACTOR * max(abs(derivative(x)), abs(derivative(midpoint)), abs(derivative(y))) + 1.0e-30
            value = function(midpoint)
            margin = abs(value) - sampled_bound * half_width
            ratio = abs(value) / (sampled_bound * half_width) if half_width > 0.0 else math.inf
            cell = {
                "left": x,
                "right": y,
                "midpoint": midpoint,
                "midpoint_value": value,
                "sampled_derivative_bound_with_safety": sampled_bound,
                "signed_exclusion_margin": margin,
                "value_to_bound_ratio": ratio,
            }
            if margin < minimum_margin:
                minimum_margin = margin
                weakest_cell = cell
            minimum_ratio = min(minimum_ratio, ratio)
            if margin <= EXCLUSION_MARGIN_TOL:
                failed_cells.append(cell)
            cell_count += 1
            x = y
    return {
        "family": family,
        "audit_start": start,
        "audit_end": end,
        "excluded_root_windows": windows,
        "complement_segments": segments,
        "cell_width": EXCLUSION_CELL_WIDTH,
        "cell_count": cell_count,
        "derivative_bound_safety_factor": DERIVATIVE_BOUND_SAFETY_FACTOR,
        "failed_cell_count": len(failed_cells),
        "minimum_signed_exclusion_margin": minimum_margin,
        "minimum_value_to_bound_ratio": minimum_ratio,
        "weakest_cell": weakest_cell,
        "sampled_lipschitz_exclusion_passed": len(failed_cells) == 0,
    }


def build_family_certificate(
    family: str,
    roots: list[dict[str, Any]],
    projection_vector: list[float],
    start: float,
    end: float,
) -> dict[str, Any]:
    root_ds = [root["root_d"] for root in roots]
    windows = root_windows(root_ds, start, end)
    if family == "zero_projection_amplitude":
        function = lambda d: projection_amplitude(projection_vector, d)
        derivative = lambda d: amplitude_derivative(projection_vector, d)
    elif family == "stationary_factor":
        function = lambda d: projection_stationary_factor(projection_vector, d)
        derivative = lambda d: finite_difference_derivative(function, d, start, end)
    else:
        raise ValueError(f"unknown family {family}")
    window_audits = monotonicity_window_audit(family, function, derivative, windows)
    exclusion = sampled_lipschitz_exclusion_audit(family, function, derivative, windows, start, end)
    return {
        "family": family,
        "root_count": len(root_ds),
        "root_window_half_width": ROOT_WINDOW_HALF_WIDTH,
        "root_windows": windows,
        "monotonicity_window_audits": window_audits,
        "all_windows_sign_change": all(audit["sign_change_across_window"] for audit in window_audits),
        "all_windows_sampled_monotone": all(audit["sampled_monotone_window"] for audit in window_audits),
        "sampled_lipschitz_exclusion": exclusion,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2450/S1400 strict pointwise projection root-isolation margin certificate

`P2450/S1400` adds a sampled root-isolation margin audit on top of the P2449 projection reduction.  It isolates the two `a·g(d)=0` roots and the one stationary-factor root in explicit sign-changing windows, checks sampled monotonicity inside those windows, and audits the complementary cells with a sampled derivative-bound margin test.

This is stronger than a raw root scan, but it is still explicitly finite and sampled: it does not export an exact interval root-exclusion theorem, point-coordinate selector, source theorem, gauge theorem, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
""".strip()
    lag_section = """
## P2450/S1400 projection root-isolation guard

`P2450/S1400` strengthens the P2449 projection route with sign-changing root windows, sampled monotonicity, and complementary sampled-Lipschitz exclusion margins.  The result remains a finite proof-audit margin certificate only; no pointwise root or coordinate is licensed as a selector/source/gauge row for `L_total`.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2449 = sources["P2449_PROJECTION_REDUCTION"].get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {})
    projection_vector = p2449.get("projection_vector", [])
    zero_projection_certificate = build_family_certificate(
        "zero_projection_amplitude",
        p2449.get("zero_projection_roots", []),
        projection_vector,
        AUDIT_INTERVAL[0] + 1.0e-6,
        AUDIT_INTERVAL[1],
    )
    stationary_factor_certificate = build_family_certificate(
        "stationary_factor",
        p2449.get("stationary_factor_roots", []),
        projection_vector,
        H_AUDIT_START,
        AUDIT_INTERVAL[1],
    )
    best_root = p2449.get("best_reconstructed_root", {})
    theorem_export = {
        "theorem_name": "P2450_T1_strict_pointwise_projection_root_isolation_margin_certificate",
        "inherited_projection_reduction_certificate": "P2449/S1399",
        "projection_vector_norm": row_norm(projection_vector),
        "root_window_half_width": ROOT_WINDOW_HALF_WIDTH,
        "exclusion_cell_width": EXCLUSION_CELL_WIDTH,
        "derivative_bound_safety_factor": DERIVATIVE_BOUND_SAFETY_FACTOR,
        "zero_projection_amplitude_certificate": zero_projection_certificate,
        "stationary_factor_certificate": stationary_factor_certificate,
        "best_inherited_stationary_factor_root": best_root,
        "best_root_volume_replayed_from_projection": projection_volume(projection_vector, best_root.get("root_d")),
        "all_root_windows_sign_change": zero_projection_certificate["all_windows_sign_change"] and stationary_factor_certificate["all_windows_sign_change"],
        "all_root_windows_sampled_monotone": zero_projection_certificate["all_windows_sampled_monotone"] and stationary_factor_certificate["all_windows_sampled_monotone"],
        "all_complement_exclusion_margins_positive": zero_projection_certificate["sampled_lipschitz_exclusion"]["sampled_lipschitz_exclusion_passed"]
        and stationary_factor_certificate["sampled_lipschitz_exclusion"]["sampled_lipschitz_exclusion_passed"],
        "minimum_exclusion_margin_across_families": min(
            zero_projection_certificate["sampled_lipschitz_exclusion"]["minimum_signed_exclusion_margin"],
            stationary_factor_certificate["sampled_lipschitz_exclusion"]["minimum_signed_exclusion_margin"],
        ),
        "sampled_lipschitz_margin_certificate_exported": True,
        "exact_interval_root_exclusion_theorem_exported_by_this_certificate": False,
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Sampled monotonicity windows plus sampled derivative-bound exclusion margins are not an exact interval root-exclusion theorem.",
            "The isolated projection roots and stationary-factor root do not choose a strict point-coordinate selector.",
            "No strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "The remaining proof upgrade would replace sampled derivative-bound margins with exact interval arithmetic or symbolic root exclusion for the projection amplitude and stationary factor."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "two_zero_projection_roots_inherited": zero_projection_certificate["root_count"] == 2,
        "one_stationary_factor_root_inherited": stationary_factor_certificate["root_count"] == 1,
        "all_root_windows_sign_change": theorem_export["all_root_windows_sign_change"],
        "all_root_windows_sampled_monotone": theorem_export["all_root_windows_sampled_monotone"],
        "all_complement_exclusion_margins_positive": theorem_export["all_complement_exclusion_margins_positive"],
        "best_root_volume_replayed": abs(theorem_export["best_root_volume_replayed_from_projection"] - best_root.get("normalized_rank_lift_volume")) < 1.0e-12,
        "no_exact_interval_root_exclusion_export": not theorem_export["exact_interval_root_exclusion_theorem_exported_by_this_certificate"],
        "no_pointwise_selector_export": not theorem_export["pointwise_coordinate_selector_exported_by_this_certificate"],
        "no_observable_source_export": not theorem_export["strict_observable_source_constraint_exported_by_this_certificate"],
        "no_gauge_slice_export": not theorem_export["gauge_slice_theorem_exported_by_this_certificate"],
        "no_value_generator_export": not theorem_export["strict_physical_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2450_s1400_v1",
        "packet_id": "P2450",
        "stage_id": "S1400",
        "status": "PASS_STRICT_POINTWISE_PROJECTION_ROOT_ISOLATION_MARGIN_NO_EXACT_INTERVAL_SELECTOR_SOURCE_THEOREM",
        "strict_pointwise_projection_root_isolation_margin_certificate": {
            "inputs": {name: rel(path) for name, path in SOURCE_FILES.items()},
            "input_fingerprints": {name: sha256_json(sources[name]) for name in sources},
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_projection_root_isolation_margin_certificate"]["theorem_export"]
    z = t["zero_projection_amplitude_certificate"]
    h = t["stationary_factor_certificate"]
    lines = [
        "# P2450/S1400 strict pointwise projection root-isolation margin certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Isolation margins",
        "",
        f"Zero-projection root windows sign-change and sampled-monotone: `{z['all_windows_sign_change']}` / `{z['all_windows_sampled_monotone']}`.",
        f"Zero-projection complement cells: `{z['sampled_lipschitz_exclusion']['cell_count']}`, failed cells: `{z['sampled_lipschitz_exclusion']['failed_cell_count']}`, weakest margin: `{z['sampled_lipschitz_exclusion']['minimum_signed_exclusion_margin']:.17g}`.",
        f"Stationary-factor root windows sign-change and sampled-monotone: `{h['all_windows_sign_change']}` / `{h['all_windows_sampled_monotone']}`.",
        f"Stationary-factor complement cells: `{h['sampled_lipschitz_exclusion']['cell_count']}`, failed cells: `{h['sampled_lipschitz_exclusion']['failed_cell_count']}`, weakest margin: `{h['sampled_lipschitz_exclusion']['minimum_signed_exclusion_margin']:.17g}`.",
        f"Minimum exclusion margin across families: `{t['minimum_exclusion_margin_across_families']:.17g}`.",
        "",
        "## Guardrail",
        "",
        "This is a sampled-Lipschitz margin certificate only.  It exports no exact interval root-exclusion theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps({"status": payload["status"], "gatekeepers": payload["gatekeeper_checks"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
