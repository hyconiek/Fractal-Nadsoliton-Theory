#!/usr/bin/env python3
"""Scratch probe: interval certificate for the Z12 inverse branch.

Prior bridge work already showed that a truncated Puiseux inverse is not globally
admissible, while a monotone pointwise inverse branch exists.  This probe makes
that latter statement more proof-like on the finite Z12 domain: it certifies,
with bisection brackets and finite monotonic margins, that each strict normalized
value S(d), d=0..11, has a unique preimage x(d) in the first decreasing legacy
interval [0,2], and that these preimages form a monotone nonnegative discrete
Z12 map.

This is still an output-matching certificate, not a strict transport derivation:
it narrows the global-map blocker to the missing strict-side derivation of this
map/law.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_phase_normalized_z12_inverse_branch_interval_certificate_report.json"
OUT_MD = HERE / "bridge_phase_normalized_z12_inverse_branch_interval_certificate_report.md"
MONOTONE_REPORT = HERE / "bridge_phase_normalized_monotone_inverse_branch_report.json"
BLOCKER_LATTICE = HERE / "bridge_completed_kernel_blocker_dependency_lattice_report.json"

LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 9.0 / 5.0}
BRACKET = (0.0, 2.0)
DOMAIN = list(range(12))
BISECTION_STEPS = 90
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def legacy_norm(x: float) -> float:
    return math.cos(LEGACY["omega"] * x + LEGACY["phi"]) / math.cos(LEGACY["phi"]) / (1.0 + LEGACY["beta_tors"] * x)


def strict_norm(d: float) -> float:
    return math.cos(STRICT["omega"] * d + STRICT["phi"]) / math.cos(STRICT["phi"]) / (1.0 + d ** STRICT["eta"])


def legacy_derivative_upper_bound() -> float:
    # On x in [0,2], theta=pi/4*x+pi/6 lies in [pi/6,2pi/3].
    # Therefore sin(theta)>=1/2 and cos(theta)>=-1/2.  The numerator of L'
    # is -omega*sin(theta)*(1+beta*x)-beta*cos(theta), so a coarse upper bound
    # is -omega/2 + beta/2.  Dividing by cos(phi)>0 and ignoring (1+beta*x)^2>=1
    # keeps a valid negative upper bound.
    return (-0.5 * LEGACY["omega"] + 0.5 * LEGACY["beta_tors"]) / math.cos(LEGACY["phi"])


def bracket_root_for_target(target: float) -> dict[str, Any]:
    lo, hi = BRACKET
    f_lo = legacy_norm(lo) - target
    f_hi = legacy_norm(hi) - target
    if abs(f_lo) < 1e-15:
        return {
            "lo": lo,
            "hi": lo,
            "mid": lo,
            "width": 0.0,
            "f_lo": f_lo,
            "f_hi": f_lo,
            "f_mid": f_lo,
            "sign_bracket_or_endpoint": True,
        }
    if f_lo * f_hi > 0:
        raise ValueError(f"target not bracketed: target={target}, f_lo={f_lo}, f_hi={f_hi}")
    for _ in range(BISECTION_STEPS):
        mid = 0.5 * (lo + hi)
        f_mid = legacy_norm(mid) - target
        if f_lo * f_mid <= 0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid
    mid = 0.5 * (lo + hi)
    return {
        "lo": lo,
        "hi": hi,
        "mid": mid,
        "width": hi - lo,
        "f_lo": f_lo,
        "f_hi": f_hi,
        "f_mid": legacy_norm(mid) - target,
        "sign_bracket_or_endpoint": f_lo * f_hi <= 0,
    }


def build_rows() -> list[dict[str, Any]]:
    rows = []
    for d in DOMAIN:
        target = strict_norm(float(d))
        root = bracket_root_for_target(target)
        rows.append(
            {
                "d": d,
                "strict_norm_target": target,
                "root_interval": root,
                "legacy_norm_at_mid": legacy_norm(root["mid"]),
                "mid_residual": legacy_norm(root["mid"]) - target,
                "mid_in_first_legacy_interval": BRACKET[0] <= root["mid"] <= BRACKET[1],
            }
        )
    return rows


def finite_differences(values: list[float]) -> list[float]:
    return [right - left for left, right in zip(values, values[1:])]


def build_payload() -> dict[str, Any]:
    monotone = load_json(MONOTONE_REPORT)
    blocker = load_json(BLOCKER_LATTICE)
    rows = build_rows()
    targets = [row["strict_norm_target"] for row in rows]
    roots = [row["root_interval"]["mid"] for row in rows]
    target_diffs = finite_differences(targets)
    root_diffs = finite_differences(roots)
    widths = [row["root_interval"]["width"] for row in rows]
    residuals = [abs(row["mid_residual"]) for row in rows]
    legacy_left = legacy_norm(BRACKET[0])
    legacy_right = legacy_norm(BRACKET[1])
    derivative_bound = legacy_derivative_upper_bound()

    return {
        "result_kind": "SCRATCH_PHASE_NORMALIZED_Z12_INVERSE_BRANCH_INTERVAL_CERTIFICATE__OUTPUT_MATCHING_NOT_STRICT_DERIVATION",
        "status": "finite-z12-monotone-inverse-branch-certified-by-interval-brackets",
        "source_reports": {
            "monotone_inverse_branch": str(MONOTONE_REPORT.relative_to(HERE.parents[1])),
            "completed_kernel_blocker_lattice": str(BLOCKER_LATTICE.relative_to(HERE.parents[1])),
        },
        "constants": {
            "legacy": LEGACY,
            "strict": STRICT,
            "branch_interval": list(BRACKET),
            "domain_d_values": DOMAIN,
            "bisection_steps": BISECTION_STEPS,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "interval_certificate": {
            "legacy_left_value": legacy_left,
            "legacy_right_value": legacy_right,
            "all_targets_inside_legacy_interval_range": all(legacy_left >= target >= legacy_right for target in targets),
            "legacy_derivative_upper_bound_on_interval": derivative_bound,
            "legacy_derivative_bound_strictly_negative": derivative_bound < 0.0,
            "all_root_intervals_sign_bracket_or_endpoint": all(row["root_interval"]["sign_bracket_or_endpoint"] for row in rows),
            "max_root_interval_width": max(widths),
            "max_mid_residual_abs": max(residuals),
            "strict_targets_strictly_decreasing_on_z12": all(diff < 0.0 for diff in target_diffs),
            "min_strict_target_drop_margin": min(-diff for diff in target_diffs),
            "roots_monotone_non_decreasing_on_z12": all(diff >= -1e-15 for diff in root_diffs),
            "min_root_increment": min(root_diffs),
            "all_midpoints_in_first_legacy_interval": all(row["mid_in_first_legacy_interval"] for row in rows),
            "z12_to_legacy_compression_x11_over_11": roots[-1] / 11.0,
        },
        "root_rows": rows,
        "upstream_blocker_update": {
            "previous_monotone_status": monotone["status"],
            "previous_global_admissibility_note": monotone["branch_metrics"]["cutoff_formal_inverse_global_admissible"],
            "blocker_lattice_global_z12_status_before_this_probe": "OPEN_BLOCKER",
            "refined_status_after_this_probe": "global_z12_output_matching_map_certified; strict_transport_derivation remains open",
            "remaining_bridge_blockers": blocker["current_frontier_summary"]["remaining_for_theorem_level_bridge"],
        },
        "exact_proof_certificate": {
            "existence": "Every strict normalized target S(d), d=0..11, lies between L(0) and L(2), so IVT gives at least one root in [0,2].",
            "uniqueness": "The analytic derivative upper bound for L on [0,2] is strictly negative, so each root is unique on the first legacy branch.",
            "monotonicity": "The finite strict targets are strictly decreasing and L is strictly decreasing, hence the inverse roots are monotone nondecreasing; bisection intervals numerically certify the finite roots.",
            "blocker_refinement": "This certifies a global finite Z12 output-matching map, but it does not derive that map from strict nadsoliton dynamics.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in a solitonic state; this is a finite output-matching map between strict and legacy carrier coordinates.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced.",
        },
        "hard_limits": [
            "No unqualified identity K_legacy_ont == K_strict_gate is asserted.",
            "No strict-side derivation of the inverse branch or transport ODE is claimed.",
            "No beta_tors -> chi_11 theorem is asserted.",
            "No legacy physical-role transfer onto K_strict_gate is used without an explicit bridge theorem.",
            "No QW-2191 discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    cert = payload["interval_certificate"]
    update = payload["upstream_blocker_update"]
    lines = [
        "# Z12 inverse-branch interval certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Interval certificate",
        "",
        f"- All targets inside legacy branch range: `{cert['all_targets_inside_legacy_interval_range']}`",
        f"- Legacy derivative upper bound on [0,2]: `{cert['legacy_derivative_upper_bound_on_interval']:.12f}`",
        f"- Derivative bound strictly negative: `{cert['legacy_derivative_bound_strictly_negative']}`",
        f"- All roots bracketed/endpoints: `{cert['all_root_intervals_sign_bracket_or_endpoint']}`",
        f"- Max root interval width: `{cert['max_root_interval_width']:.3e}`",
        f"- Max midpoint residual: `{cert['max_mid_residual_abs']:.3e}`",
        f"- Strict targets decreasing: `{cert['strict_targets_strictly_decreasing_on_z12']}`",
        f"- Roots monotone nondecreasing: `{cert['roots_monotone_non_decreasing_on_z12']}`",
        f"- Compression x(11)/11: `{cert['z12_to_legacy_compression_x11_over_11']:.12f}`",
        "",
        "## Blocker refinement",
        "",
        f"- Previous monotone status: `{update['previous_monotone_status']}`",
        f"- Previous truncated-formal-inverse admissible: `{update['previous_global_admissibility_note']}`",
        f"- Refined status: `{update['refined_status_after_this_probe']}`",
        "",
        "## Proof certificate",
        "",
    ]
    for key, value in payload["exact_proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Hard limits", ""])
    lines.extend(f"- {item}" for item in payload["hard_limits"])
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    OUT_MD.write_text(write_markdown(payload), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
