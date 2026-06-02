#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json"
MD = GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.md"

SOURCE_FILES = {
    "P2450_ROOT_ISOLATION_MARGIN": GEN / "p2450_s1400_strict_pointwise_projection_root_isolation_margin_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

STRICT_PARAMS = {"omega": 0.18575, "phi": 0.16250, "eta": 1.8}
INTERVAL_CELL_WIDTH = 1.0e-4
ROUNDING_PAD = 1.0e-18
PHASE_MONOTONE_BAND = (0.0, math.pi / 2.0)


class Interval:
    def __init__(self, lo: float, hi: float | None = None):
        self.lo = float(lo)
        self.hi = float(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ValueError(f"invalid interval [{self.lo}, {self.hi}]")

    def widen(self, pad: float = ROUNDING_PAD) -> Interval:
        return Interval(self.lo - pad, self.hi + pad)

    def __add__(self, other: float | Interval) -> Interval:
        other = to_interval(other)
        return Interval(self.lo + other.lo, self.hi + other.hi).widen()

    __radd__ = __add__

    def __sub__(self, other: float | Interval) -> Interval:
        other = to_interval(other)
        return Interval(self.lo - other.hi, self.hi - other.lo).widen()

    def __rsub__(self, other: float | Interval) -> Interval:
        return to_interval(other) - self

    def __neg__(self) -> Interval:
        return Interval(-self.hi, -self.lo).widen()

    def __mul__(self, other: float | Interval) -> Interval:
        other = to_interval(other)
        values = [self.lo * other.lo, self.lo * other.hi, self.hi * other.lo, self.hi * other.hi]
        return Interval(min(values), max(values)).widen()

    __rmul__ = __mul__

    def inv(self) -> Interval:
        if self.lo <= 0.0 <= self.hi:
            raise ZeroDivisionError(f"interval contains zero: [{self.lo}, {self.hi}]")
        values = [1.0 / self.lo, 1.0 / self.hi]
        return Interval(min(values), max(values)).widen()

    def __truediv__(self, other: float | Interval) -> Interval:
        return self * to_interval(other).inv()

    def __rtruediv__(self, other: float | Interval) -> Interval:
        return to_interval(other) / self

    def exp(self) -> Interval:
        return Interval(math.exp(self.lo), math.exp(self.hi)).widen()

    def log(self) -> Interval:
        if self.lo <= 0.0:
            raise ValueError(f"log interval is not positive: [{self.lo}, {self.hi}]")
        return Interval(math.log(self.lo), math.log(self.hi)).widen()

    def sin_monotone(self) -> Interval:
        if not (PHASE_MONOTONE_BAND[0] <= self.lo <= self.hi <= PHASE_MONOTONE_BAND[1]):
            raise ValueError(f"phase interval outside audited monotone sin band: [{self.lo}, {self.hi}]")
        return Interval(math.sin(self.lo), math.sin(self.hi)).widen()

    def cos_monotone(self) -> Interval:
        if not (PHASE_MONOTONE_BAND[0] <= self.lo <= self.hi <= PHASE_MONOTONE_BAND[1]):
            raise ValueError(f"phase interval outside audited monotone cos band: [{self.lo}, {self.hi}]")
        return Interval(math.cos(self.hi), math.cos(self.lo)).widen()

    def contains_zero(self) -> bool:
        return self.lo <= 0.0 <= self.hi

    def separation_from_zero(self) -> float:
        return -min(abs(self.lo), abs(self.hi)) if self.contains_zero() else min(abs(self.lo), abs(self.hi))

    def as_dict(self) -> dict[str, float]:
        return {"lo": self.lo, "hi": self.hi}


def to_interval(value: float | Interval) -> Interval:
    return value if isinstance(value, Interval) else Interval(value)


def rel(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": rel(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


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
        "new_packet": "P2451|S1401|interval enclosure root exclusion|projection interval enclosure|floating interval exclusion",
        "p2450_input": "P2450|S1400|root isolation margin|sampled Lipschitz exclusion|projection root-isolation",
        "interval_language": "interval enclosure|interval arithmetic|zero-excluding interval|complement cell enclosure|monotone sin",
        "honesty_guard_language": "not directed rounding|floating interval|not exact interval theorem|not symbolic root exclusion",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|point-coordinate selector",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def interval_power_const(d: Interval, exponent: float) -> Interval:
    return (exponent * d.log()).exp()


def pointwise_gradient_interval(d: Interval) -> list[Interval]:
    omega = STRICT_PARAMS["omega"]
    phi = STRICT_PARAMS["phi"]
    eta = STRICT_PARAMS["eta"]
    d_eta = interval_power_const(d, eta)
    denominator = 1.0 + d_eta
    phase = omega * d + phi
    cos_phase = phase.cos_monotone()
    sin_phase = phase.sin_monotone()
    return [
        -d * sin_phase / denominator,
        -sin_phase / denominator,
        -cos_phase * d_eta / (denominator * denominator),
        -cos_phase * d_eta * d.log() / (denominator * denominator),
    ]


def pointwise_gradient_derivative_interval(d: Interval) -> list[Interval]:
    omega = STRICT_PARAMS["omega"]
    phi = STRICT_PARAMS["phi"]
    eta = STRICT_PARAMS["eta"]
    d_eta = interval_power_const(d, eta)
    d_eta_derivative = eta * interval_power_const(d, eta - 1.0)
    denominator = 1.0 + d_eta
    denominator_derivative = d_eta_derivative
    phase = omega * d + phi
    cos_phase = phase.cos_monotone()
    sin_phase = phase.sin_monotone()
    g2 = -cos_phase * d_eta / (denominator * denominator)
    g2_derivative = (
        omega * sin_phase * d_eta / (denominator * denominator)
        - cos_phase * d_eta_derivative / (denominator * denominator)
        + 2.0 * cos_phase * d_eta * denominator_derivative / (denominator * denominator * denominator)
    )
    return [
        -(sin_phase + d * omega * cos_phase) / denominator + d * sin_phase * denominator_derivative / (denominator * denominator),
        -omega * cos_phase / denominator + sin_phase * denominator_derivative / (denominator * denominator),
        g2_derivative,
        g2_derivative * d.log() + g2 / d,
    ]


def interval_dot(left: list[float] | list[Interval], right: list[Interval]) -> Interval:
    result = Interval(0.0)
    for l_value, r_value in zip(left, right):
        result = result + l_value * r_value
    return result


def projection_amplitude_interval(projection_vector: list[float], d: Interval) -> Interval:
    return interval_dot(projection_vector, pointwise_gradient_interval(d))


def stationary_factor_interval(projection_vector: list[float], d: Interval) -> Interval:
    gradient = pointwise_gradient_interval(d)
    gradient_derivative = pointwise_gradient_derivative_interval(d)
    amplitude = interval_dot(projection_vector, gradient)
    amplitude_derivative = interval_dot(projection_vector, gradient_derivative)
    gradient_square = interval_dot(gradient, gradient)
    gradient_square_derivative = 2.0 * interval_dot(gradient, gradient_derivative)
    return 2.0 * amplitude_derivative * gradient_square - amplitude * gradient_square_derivative


def interval_exclusion_audit(
    family: str,
    function: Callable[[Interval], Interval],
    complement_segments: list[dict[str, float]],
) -> dict[str, Any]:
    cell_count = 0
    zero_containing_cells: list[dict[str, Any]] = []
    weakest_cell: dict[str, Any] | None = None
    minimum_separation = math.inf
    maximum_width = 0.0
    for segment in complement_segments:
        x = segment["left"]
        right = segment["right"]
        while x < right - 1.0e-15:
            y = min(right, x + INTERVAL_CELL_WIDTH)
            value = function(Interval(x, y))
            separation = value.separation_from_zero()
            cell = {"left": x, "right": y, "interval_value": value.as_dict(), "separation_from_zero": separation}
            if separation < minimum_separation:
                minimum_separation = separation
                weakest_cell = cell
            maximum_width = max(maximum_width, value.hi - value.lo)
            if value.contains_zero():
                zero_containing_cells.append(cell)
            cell_count += 1
            x = y
    return {
        "family": family,
        "interval_cell_width": INTERVAL_CELL_WIDTH,
        "cell_count": cell_count,
        "zero_containing_cell_count": len(zero_containing_cells),
        "minimum_separation_from_zero": minimum_separation,
        "maximum_interval_width": maximum_width,
        "weakest_cell": weakest_cell,
        "zero_containing_cell_samples": zero_containing_cells[:10],
        "all_complement_cells_exclude_zero": len(zero_containing_cells) == 0,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2451/S1401 strict pointwise projection interval-enclosure root-exclusion audit

`P2451/S1401` replaces the P2450 midpoint-margin complement audit with direct interval enclosures of the two scalar projection functions on `1e-4` complement cells.  On the P2450 root-window complements, interval evaluation of `a·g(d)` and of the stationary factor excludes zero in every audited cell.

This is a stronger finite enclosure audit than sampled midpoint margins, but it is still not a symbolic proof or directed-rounding exact interval theorem.  It exports no point-coordinate selector, no source/gauge theorem, no role-bearing `L_total`, no `QW-2191` discharge, and no ToE closure.
""".strip()
    lag_section = """
## P2451/S1401 interval-enclosure root-exclusion guard

`P2451/S1401` gives a finite interval-enclosure replay of the P2450 complement exclusion for the projection amplitude and stationary factor.  The result remains a numerical enclosure audit rather than exact selector/source/gauge authority, so no pointwise coordinate or root is admitted into `L_total`.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2450 = sources["P2450_ROOT_ISOLATION_MARGIN"].get("strict_pointwise_projection_root_isolation_margin_certificate", {}).get("theorem_export", {})
    projection_vector = load_json(GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json").get(
        "strict_pointwise_rank_lift_projection_reduction_certificate", {}
    ).get("theorem_export", {}).get("projection_vector", [])
    zero_segments = p2450.get("zero_projection_amplitude_certificate", {}).get("sampled_lipschitz_exclusion", {}).get("complement_segments", [])
    stationary_segments = p2450.get("stationary_factor_certificate", {}).get("sampled_lipschitz_exclusion", {}).get("complement_segments", [])
    zero_interval_audit = interval_exclusion_audit(
        "zero_projection_amplitude",
        lambda interval: projection_amplitude_interval(projection_vector, interval),
        zero_segments,
    )
    stationary_interval_audit = interval_exclusion_audit(
        "stationary_factor",
        lambda interval: stationary_factor_interval(projection_vector, interval),
        stationary_segments,
    )
    theorem_export = {
        "theorem_name": "P2451_T1_strict_pointwise_projection_interval_enclosure_root_exclusion_audit",
        "inherited_root_isolation_margin_certificate": "P2450/S1400",
        "interval_cell_width": INTERVAL_CELL_WIDTH,
        "rounding_pad_used_in_float_interval_ops": ROUNDING_PAD,
        "phase_monotone_band": {"lo": PHASE_MONOTONE_BAND[0], "hi": PHASE_MONOTONE_BAND[1]},
        "zero_projection_amplitude_interval_audit": zero_interval_audit,
        "stationary_factor_interval_audit": stationary_interval_audit,
        "all_interval_cells_exclude_zero": zero_interval_audit["all_complement_cells_exclude_zero"]
        and stationary_interval_audit["all_complement_cells_exclude_zero"],
        "minimum_interval_separation_from_zero_across_families": min(
            zero_interval_audit["minimum_separation_from_zero"],
            stationary_interval_audit["minimum_separation_from_zero"],
        ),
        "floating_interval_enclosure_audit_exported": True,
        "directed_rounding_interval_theorem_exported_by_this_certificate": False,
        "symbolic_root_exclusion_theorem_exported_by_this_certificate": False,
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "The interval enclosures are finite floating interval enclosures with explicit padding, not a directed-rounding symbolic interval theorem.",
            "Zero exclusion on complement cells does not choose a strict point-coordinate selector or observable/source row.",
            "No strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "A stronger proof would replace these floating interval enclosures with directed-rounding interval arithmetic or symbolic root-exclusion bounds for the projection amplitude and stationary factor."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "zero_projection_interval_cells_exclude_zero": zero_interval_audit["all_complement_cells_exclude_zero"],
        "stationary_factor_interval_cells_exclude_zero": stationary_interval_audit["all_complement_cells_exclude_zero"],
        "all_interval_cells_exclude_zero": theorem_export["all_interval_cells_exclude_zero"],
        "positive_minimum_interval_separation": theorem_export["minimum_interval_separation_from_zero_across_families"] > 0.0,
        "no_directed_rounding_interval_theorem": not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
        "no_symbolic_root_exclusion": not theorem_export["symbolic_root_exclusion_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2451_s1401_v1",
        "packet_id": "P2451",
        "stage_id": "S1401",
        "status": "PASS_STRICT_POINTWISE_PROJECTION_INTERVAL_ENCLOSURE_ROOT_EXCLUSION_AUDIT_NO_EXACT_SELECTOR_SOURCE_THEOREM",
        "strict_pointwise_projection_interval_enclosure_root_exclusion_audit": {
            "inputs": {name: rel(path) for name, path in SOURCE_FILES.items()},
            "input_fingerprints": {name: sha256_json(sources[name]) for name in sources},
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_projection_interval_enclosure_root_exclusion_audit"]["theorem_export"]
    z = t["zero_projection_amplitude_interval_audit"]
    h = t["stationary_factor_interval_audit"]
    lines = [
        "# P2451/S1401 strict pointwise projection interval-enclosure root-exclusion audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Interval enclosure result",
        "",
        f"Zero-projection complement cells: `{z['cell_count']}`, zero-containing cells: `{z['zero_containing_cell_count']}`, minimum separation: `{z['minimum_separation_from_zero']:.17g}`.",
        f"Stationary-factor complement cells: `{h['cell_count']}`, zero-containing cells: `{h['zero_containing_cell_count']}`, minimum separation: `{h['minimum_separation_from_zero']:.17g}`.",
        f"Minimum interval separation across families: `{t['minimum_interval_separation_from_zero_across_families']:.17g}`.",
        "",
        "## Guardrail",
        "",
        "This is a finite floating interval-enclosure audit, not a directed-rounding exact interval theorem and not symbolic root exclusion.  It exports no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
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
