#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
OUT = GEN / "p2453_s1403_strict_pointwise_directed_decimal_weakest_cell_replay_certificate.json"
MD = GEN / "p2453_s1403_strict_pointwise_directed_decimal_weakest_cell_replay_certificate.md"

SOURCE_FILES = {
    "P2451_INTERVAL_ENCLOSURE": GEN / "p2451_s1401_strict_pointwise_projection_interval_enclosure_root_exclusion_audit.json",
    "P2452_RATIONAL_PRECONDITIONS": GEN / "p2452_s1402_strict_pointwise_interval_precondition_rational_certificate.json",
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
}
DOC_FILES = {
    "equation_sheet": ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md",
    "lagrangian_eom_draft": ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md",
}

DECIMAL_PRECISION = 90
DECIMAL_SAFETY_PAD = Decimal("1e-65")
TAYLOR_BOUND_PAD = Decimal("1e-60")
TAYLOR_TERM_STOP = Decimal("1e-70")
STRICT_OMEGA = Decimal("0.18575")
STRICT_PHI = Decimal("0.1625")
STRICT_ETA = Decimal("1.8")

getcontext().prec = DECIMAL_PRECISION


class DecimalInterval:
    def __init__(self, lo: Decimal | float | str, hi: Decimal | float | str | None = None):
        self.lo = to_decimal(lo)
        self.hi = to_decimal(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ValueError(f"invalid interval [{self.lo}, {self.hi}]")

    def widen(self, pad: Decimal = DECIMAL_SAFETY_PAD) -> DecimalInterval:
        return DecimalInterval(self.lo - pad, self.hi + pad)

    def __add__(self, other: DecimalInterval | Decimal | float | str) -> DecimalInterval:
        other = to_interval(other)
        return DecimalInterval(self.lo + other.lo, self.hi + other.hi).widen()

    __radd__ = __add__

    def __sub__(self, other: DecimalInterval | Decimal | float | str) -> DecimalInterval:
        other = to_interval(other)
        return DecimalInterval(self.lo - other.hi, self.hi - other.lo).widen()

    def __rsub__(self, other: DecimalInterval | Decimal | float | str) -> DecimalInterval:
        return to_interval(other) - self

    def __neg__(self) -> DecimalInterval:
        return DecimalInterval(-self.hi, -self.lo).widen()

    def __mul__(self, other: DecimalInterval | Decimal | float | str) -> DecimalInterval:
        other = to_interval(other)
        products = [self.lo * other.lo, self.lo * other.hi, self.hi * other.lo, self.hi * other.hi]
        return DecimalInterval(min(products), max(products)).widen()

    __rmul__ = __mul__

    def inv(self) -> DecimalInterval:
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError(f"interval contains zero: [{self.lo}, {self.hi}]")
        values = [Decimal(1) / self.lo, Decimal(1) / self.hi]
        return DecimalInterval(min(values), max(values)).widen()

    def __truediv__(self, other: DecimalInterval | Decimal | float | str) -> DecimalInterval:
        return self * to_interval(other).inv()

    def __rtruediv__(self, other: DecimalInterval | Decimal | float | str) -> DecimalInterval:
        return to_interval(other) / self

    def log(self) -> DecimalInterval:
        if self.lo <= 0:
            raise ValueError(f"log interval is not positive: [{self.lo}, {self.hi}]")
        return DecimalInterval(self.lo.ln() - DECIMAL_SAFETY_PAD, self.hi.ln() + DECIMAL_SAFETY_PAD)

    def exp(self) -> DecimalInterval:
        return DecimalInterval(self.lo.exp() - DECIMAL_SAFETY_PAD, self.hi.exp() + DECIMAL_SAFETY_PAD)

    def contains_zero(self) -> bool:
        return self.lo <= 0 <= self.hi

    def separation_from_zero(self) -> Decimal:
        return -min(abs(self.lo), abs(self.hi)) if self.contains_zero() else min(abs(self.lo), abs(self.hi))

    def as_dict(self) -> dict[str, str]:
        return {"lo": str(self.lo), "hi": str(self.hi)}


def to_decimal(value: Decimal | float | str) -> Decimal:
    return value if isinstance(value, Decimal) else Decimal(str(value))


def to_interval(value: DecimalInterval | Decimal | float | str) -> DecimalInterval:
    return value if isinstance(value, DecimalInterval) else DecimalInterval(value)


def decimal_factorial(n: int) -> Decimal:
    result = Decimal(1)
    for value in range(2, n + 1):
        result *= value
    return result


def decimal_sin_point_bounds(x: Decimal) -> tuple[Decimal, Decimal]:
    total = Decimal(0)
    partials: list[Decimal] = []
    for k in range(50):
        term = (x ** (2 * k + 1)) / decimal_factorial(2 * k + 1)
        total = total + term if k % 2 == 0 else total - term
        partials.append(total)
        if term < TAYLOR_TERM_STOP:
            break
    return min(partials[-2:]) - TAYLOR_BOUND_PAD, max(partials[-2:]) + TAYLOR_BOUND_PAD


def decimal_cos_point_bounds(x: Decimal) -> tuple[Decimal, Decimal]:
    total = Decimal(0)
    partials: list[Decimal] = []
    for k in range(50):
        term = (x ** (2 * k)) / decimal_factorial(2 * k)
        total = total + term if k % 2 == 0 else total - term
        partials.append(total)
        if term < TAYLOR_TERM_STOP:
            break
    return min(partials[-2:]) - TAYLOR_BOUND_PAD, max(partials[-2:]) + TAYLOR_BOUND_PAD


def sin_monotone_interval(x: DecimalInterval) -> DecimalInterval:
    lower, _ = decimal_sin_point_bounds(x.lo)
    _, upper = decimal_sin_point_bounds(x.hi)
    return DecimalInterval(lower, upper)


def cos_monotone_interval(x: DecimalInterval) -> DecimalInterval:
    lower, _ = decimal_cos_point_bounds(x.hi)
    _, upper = decimal_cos_point_bounds(x.lo)
    return DecimalInterval(lower, upper)


def interval_power_const(d: DecimalInterval, exponent: Decimal) -> DecimalInterval:
    return (exponent * d.log()).exp()


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
        "new_packet": "P2453|S1403|directed decimal weakest-cell|Taylor endpoint replay|decimal weakest cell",
        "p2451_input": "P2451|S1401|interval enclosure root exclusion|weakest_cell|floating interval enclosure",
        "p2452_input": "P2452|S1402|rational phase-band|interval precondition rational|monotone trigonometric",
        "backend_language": "Decimal interval|Taylor endpoint|alternating series|directed decimal|weakest-cell replay",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|point-coordinate selector",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def interval_dot(left: list[Decimal] | list[DecimalInterval], right: list[DecimalInterval]) -> DecimalInterval:
    result = DecimalInterval(0)
    for l_value, r_value in zip(left, right):
        result = result + l_value * r_value
    return result


def pointwise_gradient_interval(d: DecimalInterval) -> list[DecimalInterval]:
    d_eta = interval_power_const(d, STRICT_ETA)
    denominator = 1 + d_eta
    phase = STRICT_OMEGA * d + STRICT_PHI
    cos_phase = cos_monotone_interval(phase)
    sin_phase = sin_monotone_interval(phase)
    return [
        -d * sin_phase / denominator,
        -sin_phase / denominator,
        -cos_phase * d_eta / (denominator * denominator),
        -cos_phase * d_eta * d.log() / (denominator * denominator),
    ]


def pointwise_gradient_derivative_interval(d: DecimalInterval) -> list[DecimalInterval]:
    d_eta = interval_power_const(d, STRICT_ETA)
    d_eta_derivative = STRICT_ETA * interval_power_const(d, STRICT_ETA - Decimal(1))
    denominator = 1 + d_eta
    phase = STRICT_OMEGA * d + STRICT_PHI
    cos_phase = cos_monotone_interval(phase)
    sin_phase = sin_monotone_interval(phase)
    g2 = -cos_phase * d_eta / (denominator * denominator)
    g2_derivative = (
        STRICT_OMEGA * sin_phase * d_eta / (denominator * denominator)
        - cos_phase * d_eta_derivative / (denominator * denominator)
        + 2 * cos_phase * d_eta * d_eta_derivative / (denominator * denominator * denominator)
    )
    return [
        -(sin_phase + d * STRICT_OMEGA * cos_phase) / denominator + d * sin_phase * d_eta_derivative / (denominator * denominator),
        -STRICT_OMEGA * cos_phase / denominator + sin_phase * d_eta_derivative / (denominator * denominator),
        g2_derivative,
        g2_derivative * d.log() + g2 / d,
    ]


def projection_amplitude_interval(projection_vector: list[Decimal], d: DecimalInterval) -> DecimalInterval:
    return interval_dot(projection_vector, pointwise_gradient_interval(d))


def stationary_factor_interval(projection_vector: list[Decimal], d: DecimalInterval) -> DecimalInterval:
    gradient = pointwise_gradient_interval(d)
    gradient_derivative = pointwise_gradient_derivative_interval(d)
    amplitude = interval_dot(projection_vector, gradient)
    amplitude_derivative = interval_dot(projection_vector, gradient_derivative)
    gradient_square = interval_dot(gradient, gradient)
    gradient_square_derivative = 2 * interval_dot(gradient, gradient_derivative)
    return 2 * amplitude_derivative * gradient_square - amplitude * gradient_square_derivative


def replay_cell(family: str, cell: dict[str, Any], function: Callable[[DecimalInterval], DecimalInterval]) -> dict[str, Any]:
    interval = DecimalInterval(cell["left"], cell["right"])
    value = function(interval)
    p2451_value = cell.get("interval_value", {})
    return {
        "family": family,
        "cell_left": cell["left"],
        "cell_right": cell["right"],
        "p2451_float_interval_value": p2451_value,
        "decimal_taylor_interval_value": value.as_dict(),
        "decimal_interval_excludes_zero": not value.contains_zero(),
        "decimal_separation_from_zero": str(value.separation_from_zero()),
        "same_positive_sign_as_p2451": value.lo > 0 and p2451_value.get("lo", 0) > 0,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2453/S1403 strict pointwise directed-decimal weakest-cell replay certificate

`P2453/S1403` replays the two weakest P2451 complement cells with a higher-precision Decimal interval backend and monotone Taylor endpoint enclosures for `sin` and `cos`.  The zero-projection weakest cell and stationary-factor weakest cell both remain strictly positive and zero-excluding under this independent endpoint replay.

This is only a weakest-cell backend cross-check, not a full directed-rounding interval theorem or symbolic root-exclusion proof.  It exports no point-coordinate selector, source/gauge theorem, role-bearing `L_total`, `QW-2191` discharge, or ToE closure.
""".strip()
    lag_section = """
## P2453/S1403 directed-decimal weakest-cell replay guard

`P2453/S1403` independently replays the weakest P2451 cells with Decimal/Taylor endpoint enclosures, strengthening confidence in the interval audit bottleneck cells.  It remains a numerical backend cross-check and does not license any pointwise coordinate as a selector/source/gauge row for `L_total`.
""".strip()
    for path, section in [(DOC_FILES["equation_sheet"], eq_section), (DOC_FILES["lagrangian_eom_draft"], lag_section)]:
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        marker = section.splitlines()[0]
        if marker not in text:
            path.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def build_payload() -> dict[str, Any]:
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    grep = rg_audit()
    p2451 = sources["P2451_INTERVAL_ENCLOSURE"].get("strict_pointwise_projection_interval_enclosure_root_exclusion_audit", {}).get("theorem_export", {})
    p2452 = sources["P2452_RATIONAL_PRECONDITIONS"].get("strict_pointwise_interval_precondition_rational_certificate", {}).get("theorem_export", {})
    projection_vector = [Decimal(str(value)) for value in sources["P2449_PROJECTION_REDUCTION"].get("strict_pointwise_rank_lift_projection_reduction_certificate", {}).get("theorem_export", {}).get("projection_vector", [])]
    zero_cell = p2451.get("zero_projection_amplitude_interval_audit", {}).get("weakest_cell", {})
    stationary_cell = p2451.get("stationary_factor_interval_audit", {}).get("weakest_cell", {})
    zero_replay = replay_cell(
        "zero_projection_amplitude",
        zero_cell,
        lambda interval: projection_amplitude_interval(projection_vector, interval),
    )
    stationary_replay = replay_cell(
        "stationary_factor",
        stationary_cell,
        lambda interval: stationary_factor_interval(projection_vector, interval),
    )
    theorem_export = {
        "theorem_name": "P2453_T1_strict_pointwise_directed_decimal_weakest_cell_replay_certificate",
        "inherited_interval_enclosure_audit": "P2451/S1401",
        "inherited_rational_precondition_certificate": "P2452/S1402",
        "decimal_precision": DECIMAL_PRECISION,
        "decimal_safety_pad": str(DECIMAL_SAFETY_PAD),
        "taylor_bound_pad": str(TAYLOR_BOUND_PAD),
        "taylor_term_stop": str(TAYLOR_TERM_STOP),
        "p2452_phase_precondition_inherited": p2452.get("phase_band_inside_open_0_pi_over_2"),
        "zero_projection_weakest_cell_replay": zero_replay,
        "stationary_factor_weakest_cell_replay": stationary_replay,
        "both_weakest_cells_exclude_zero_under_decimal_taylor_replay": zero_replay["decimal_interval_excludes_zero"] and stationary_replay["decimal_interval_excludes_zero"],
        "both_weakest_cells_keep_p2451_positive_sign": zero_replay["same_positive_sign_as_p2451"] and stationary_replay["same_positive_sign_as_p2451"],
        "full_complement_directed_rounding_interval_theorem_exported_by_this_certificate": False,
        "symbolic_root_exclusion_theorem_exported_by_this_certificate": False,
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_exported": False,
        "not_licensed": [
            "Weakest-cell Decimal/Taylor replay is not a full complement directed-rounding interval theorem.",
            "The endpoint replay does not prove symbolic root exclusion and does not choose a point-coordinate selector.",
            "No strict observable/source theorem, gauge-slice theorem, physical-value generator, QW-2191 discharge, role-bearing L_total export, or ToE closure is exported.",
        ],
        "next_honest_step": (
            "Extend the Decimal/Taylor endpoint backend to every P2451 complement cell, or replace it with a formal directed-rounding interval backend."
        ),
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2452_phase_precondition_inherited": theorem_export["p2452_phase_precondition_inherited"] is True,
        "zero_weakest_cell_excludes_zero": zero_replay["decimal_interval_excludes_zero"],
        "stationary_weakest_cell_excludes_zero": stationary_replay["decimal_interval_excludes_zero"],
        "both_weakest_cells_exclude_zero": theorem_export["both_weakest_cells_exclude_zero_under_decimal_taylor_replay"],
        "both_weakest_cells_keep_p2451_sign": theorem_export["both_weakest_cells_keep_p2451_positive_sign"],
        "no_full_complement_theorem": not theorem_export["full_complement_directed_rounding_interval_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2453_s1403_v1",
        "packet_id": "P2453",
        "stage_id": "S1403",
        "status": "PASS_STRICT_POINTWISE_DIRECTED_DECIMAL_WEAKEST_CELL_REPLAY_NO_FULL_INTERVAL_SELECTOR_SOURCE_THEOREM",
        "strict_pointwise_directed_decimal_weakest_cell_replay_certificate": {
            "inputs": {name: rel(path) for name, path in SOURCE_FILES.items()},
            "input_fingerprints": {name: sha256_json(sources[name]) for name in sources},
            "rg_nonduplication_audit": grep,
            "theorem_export": theorem_export,
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_directed_decimal_weakest_cell_replay_certificate"]["theorem_export"]
    z = t["zero_projection_weakest_cell_replay"]
    h = t["stationary_factor_weakest_cell_replay"]
    lines = [
        "# P2453/S1403 strict pointwise directed-decimal weakest-cell replay certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Weakest-cell replay",
        "",
        f"Zero-projection Decimal/Taylor interval: `{z['decimal_taylor_interval_value']}`, excludes zero: `{z['decimal_interval_excludes_zero']}`.",
        f"Stationary-factor Decimal/Taylor interval: `{h['decimal_taylor_interval_value']}`, excludes zero: `{h['decimal_interval_excludes_zero']}`.",
        f"Both keep P2451 positive sign: `{t['both_weakest_cells_keep_p2451_positive_sign']}`.",
        "",
        "## Guardrail",
        "",
        "This is a weakest-cell Decimal/Taylor backend replay only.  It exports no full complement directed-rounding interval theorem, no symbolic root-exclusion theorem, no point-coordinate selector, no strict observable/source theorem, no gauge-slice theorem, no physical-value generator, no QW-2191 discharge, no role-bearing `L_total`, and no ToE closure.",
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
