#!/usr/bin/env python3
"""Execute FIN Research Programs P379--P393 (Release 10.33).

The round reflects compact exact certificate predicates through Lean, closes
the oscillatory primal-dual gap below 1e-7, proves global injectivity of the
strict radial law on integer distances, tests exactly one damping-completion
naturality square, classifies the explicit component-phase lattice of the
photonic mesh, constructs a lossy/dephasing adaptive-comb sandwich, propagates
clock tubes, builds one explicit Jordan-sampling realization functor, audits
resource conversions, and proves a dilation-semigroup no-go for a canonical
outer scale.

Programs requiring independent hardware, hold-outs, standards, reservoir
records, or electroweak unblinding remain external admission gates.
"""

from __future__ import annotations

import csv
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm

import fin_programs_255_266 as core
import fin_programs_295_308 as p295
import fin_programs_351_364 as p351
import fin_programs_365_378 as p365


ROOT = Path(__file__).resolve().parent
FIGURE_DIR = ROOT / "FIN_Programs_379_393_Figures"
RESULTS_PATH = ROOT / "FIN_Programs_379_393_Results.json"
SUMMARY_PATH = ROOT / "FIN_Programs_379_393_Summary.csv"
REFLECTION_PATH = ROOT / "FIN_Programs_379_393_Arithmetic_Reflection.csv"
CONTACT_PATH = ROOT / "FIN_Programs_379_393_Oscillatory_Contact.csv"
INJECTIVITY_PATH = ROOT / "FIN_Programs_379_393_Distance_Injectivity.csv"
DAMPING_PATH = ROOT / "FIN_Programs_379_393_Damping_Naturality.csv"
GAUGE_PATH = ROOT / "FIN_Programs_379_393_Component_Gauge.csv"
PHOTONIC_PATH = ROOT / "FIN_Programs_379_393_Photonic_Gate.csv"
LOSSY_COMB_PATH = ROOT / "FIN_Programs_379_393_Lossy_Comb.csv"
CLOCK_PATH = ROOT / "FIN_Programs_379_393_Clock_Tube.csv"
REALIZATION_PATH = ROOT / "FIN_Programs_379_393_Jordan_Realization.csv"
CONVERSION_PATH = ROOT / "FIN_Programs_379_393_Conversion_Matrix.csv"
OUTER_SCALE_PATH = ROOT / "FIN_Programs_379_393_Outer_Scale.csv"
HOLDOUT_PATH = ROOT / "FIN_Programs_379_393_Holdout_Gate.csv"
STANDARDS_PATH = ROOT / "FIN_Programs_379_393_Standards_Gate.csv"
RESERVOIR_PATH = ROOT / "FIN_Programs_379_393_Reservoir_Gate.csv"
EW_PATH = ROOT / "FIN_Programs_379_393_EW_Gate.csv"
FORMAL_SOURCE = ROOT / "FIN_Programs_379_393_Formal_Core.lean"
ARITHMETIC_SOURCE = ROOT / "FIN_Programs_379_393_Arithmetic_Reflection.lean"
LEAN_BINARY = (
    ROOT
    / ".elan"
    / "toolchains"
    / "leanprover--lean4---v4.28.0"
    / "bin"
    / "lean"
)

SEED = 20260731 + 33
N = 12
sys.set_int_max_str_digits(0)


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: (
                        json.dumps(json_ready(value), ensure_ascii=False)
                        if isinstance(value, (dict, list, tuple, np.ndarray))
                        else value
                    )
                    for key, value in row.items()
                }
            )


# ---------------------------------------------------------------------------
# Shared exact P366 primal reconstruction
# ---------------------------------------------------------------------------


def p366_primal_certificate() -> dict[str, Any]:
    """Recompute the exact seven-atom Krawczyk certificate from P366."""

    fixed_first = Fraction("0.019826532267005523")
    centers = [
        Fraction(value)
        for value in (
            "0.12942143496593873047268675338102721084490827345483",
            "0.29504200366168479087247799964895432592823041366796",
            "0.52668984536238593046432125615043289896573399940205",
            "0.81405904546345298302114172056788422010765063455263",
            "0.94184864799390211194051877479072060685679704018544",
            "-0.62639752027224283039953898685936809372793982479762",
            "0.61994981683483767070064446385072926501050616538185",
            "0.59211987331704356442867078651436937457184111020374",
            "0.32113508631461340294581553652363912311958171218392",
            "0.14195366099001397263367992494077790131159961634391",
            "-0.080955948127583162253620475604772753834474631945205",
            "0.019020934133646015552385944844570959558184305726502",
        )
    ]
    radius = Fraction(1, 10**20)
    box = [p351.RI(center - radius, center + radius) for center in centers]
    point = [p351.RI.point(center) for center in centers]
    moments = [p351.oscillatory_moment_interval(order) for order in range(12)]

    def system(variables: list[p351.RI]) -> list[p351.RI]:
        roots = [
            p351.RI.point(fixed_first),
            *variables[:5],
            p351.RI.point(1),
        ]
        weights = variables[5:]
        return [
            p351.interval_sum(
                [
                    weights[index] * roots[index].power(order)
                    for index in range(7)
                ]
            )
            - moments[order]
            for order in range(12)
        ]

    def jacobian(variables: list[p351.RI]) -> list[list[p351.RI]]:
        roots = [
            p351.RI.point(fixed_first),
            *variables[:5],
            p351.RI.point(1),
        ]
        weights = variables[5:]
        rows = []
        for order in range(12):
            root_columns = [
                (
                    weights[index]
                    * roots[index].power(order - 1)
                    * p351.RI.point(order)
                    if order
                    else p351.RI.point(0)
                )
                for index in range(1, 6)
            ]
            rows.append(
                root_columns
                + [roots[index].power(order) for index in range(7)]
            )
        return rows

    inverse = p351.exact_matrix_inverse(
        [[entry.mid for entry in row] for row in jacobian(point)]
    )
    correction = p351.matrix_vector_interval(inverse, system(point))
    base = [
        p351.RI.point(centers[index]) - correction[index]
        for index in range(12)
    ]
    cj = p351.interval_matrix_product(inverse, jacobian(box))
    delta = [p351.RI(-radius, radius) for _ in range(12)]
    image = []
    for row in range(12):
        remainder = p351.interval_sum(
            [
                (p351.RI.point(int(row == column)) - cj[row][column])
                * delta[column]
                for column in range(12)
            ]
        )
        image.append(base[row] + remainder)

    inside = [
        image[index].strictly_inside(box[index]) for index in range(12)
    ]
    signs = [
        1 if value.lo > 0 else -1 if value.hi < 0 else 0
        for value in image[5:]
    ]
    assert all(inside)
    assert signs == [-1, 1, 1, 1, 1, -1, 1]
    primal = -image[5] - image[10]
    return {
        "fixed_first": fixed_first,
        "centers": centers,
        "box": box,
        "image": image,
        "inside": inside,
        "signs": signs,
        "primal": primal,
        "moments": moments,
    }


def bernstein_to_power(coefficients: list[Fraction]) -> list[Fraction]:
    degree = len(coefficients) - 1
    power = [Fraction(0)] * (degree + 1)
    for index, coefficient in enumerate(coefficients):
        for offset in range(degree - index + 1):
            exponent = index + offset
            power[exponent] += (
                coefficient
                * math.comb(degree, index)
                * math.comb(degree - index, offset)
                * ((-1) ** offset)
            )
    return power


def p380_dual_certificate(primal_data: dict[str, Any]) -> dict[str, Any]:
    """Exact depth-14 Bernstein certificate from a frozen optimized basis."""

    raw_bernstein = [
        Fraction(value)
        for value in (
            "-0.6168528203337204",
            "-4.750518040543586",
            "22.349279010030138",
            "-59.63065668901719",
            "111.43808077322346",
            "-152.55414563959658",
            "151.2832745390066",
            "-103.5741356911199",
            "42.488230990350125",
            "-6.419183903793957",
            "-2.6492983100850838",
            "0",
        )
    ]
    depth = 14
    cells = [raw_bernstein]
    for _ in range(depth):
        cells = [
            child
            for cell in cells
            for child in p365.split_bernstein_half(cell)
        ]
    raw_lower = min(min(cell) for cell in cells)
    raw_upper = max(max(cell) for cell in cells)
    span = raw_upper - raw_lower
    safe_bernstein = [
        (coefficient - raw_upper) / span for coefficient in raw_bernstein
    ]
    safe_cells = [safe_bernstein]
    for _ in range(depth):
        safe_cells = [
            child
            for cell in safe_cells
            for child in p365.split_bernstein_half(cell)
        ]
    safe_lower = min(min(cell) for cell in safe_cells)
    safe_upper = max(max(cell) for cell in safe_cells)
    assert safe_lower >= -1
    assert safe_upper <= 0
    power = bernstein_to_power(safe_bernstein)
    dual = p351.interval_sum(
        [
            primal_data["moments"][order].scale(power[order])
            for order in range(12)
        ]
    )
    primal = primal_data["primal"]
    gap = primal.hi - dual.lo
    assert dual.lo <= primal.hi
    assert gap < Fraction(1, 10**7)
    return {
        "depth": depth,
        "raw_bernstein": raw_bernstein,
        "raw_lower": raw_lower,
        "raw_upper": raw_upper,
        "safe_bernstein": safe_bernstein,
        "safe_lower": safe_lower,
        "safe_upper": safe_upper,
        "power": power,
        "dual": dual,
        "primal": primal,
        "gap": gap,
    }


# ---------------------------------------------------------------------------
# P379: compact exact end-certificate reflection through Lean
# ---------------------------------------------------------------------------


def int_term(value: int) -> str:
    return f"({value} : Int)"


def comparison_expression(
    left: Fraction, operator: str, right: Fraction
) -> str:
    left_cross = f"{int_term(left.numerator)} * {int_term(right.denominator)}"
    right_cross = f"{int_term(right.numerator)} * {int_term(left.denominator)}"
    return f"{left_cross} {operator} {right_cross}"


def program_379(
    primal_data: dict[str, Any], dual_data: dict[str, Any]
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    claims: list[tuple[str, Fraction, str, Fraction]] = []
    for index, (box, image) in enumerate(
        zip(primal_data["box"], primal_data["image"])
    ):
        claims.append((f"krawczyk_lower_{index}", box.lo, "<", image.lo))
        claims.append((f"krawczyk_upper_{index}", image.hi, "<", box.hi))
    for index, interval in enumerate(primal_data["image"][5:]):
        if interval.lo > 0:
            claims.append((f"weight_positive_{index}", Fraction(0), "<", interval.lo))
        else:
            claims.append((f"weight_negative_{index}", interval.hi, "<", Fraction(0)))
    claims.extend(
        [
            ("bernstein_lower", Fraction(-1), "<=", dual_data["safe_lower"]),
            ("bernstein_upper", dual_data["safe_upper"], "<=", Fraction(0)),
            (
                "weak_duality_order",
                dual_data["dual"].lo,
                "<=",
                dual_data["primal"].hi,
            ),
            (
                "contact_gap_below_1e_7",
                dual_data["gap"],
                "<",
                Fraction(1, 10**7),
            ),
        ]
    )
    lines = [
        "import Std",
        "",
        "/- Generated compact arithmetic reflection for P379.",
        "   It checks exact terminal integer predicates, not the correctness",
        "   of the Python interval/Taylor/Bernstein generators. -/",
        "namespace FINPrograms379To393Arithmetic",
        "",
    ]
    rows = []
    for index, (name, left, operator, right) in enumerate(claims):
        expression = comparison_expression(left, operator, right)
        lines.extend(
            [
                f"theorem claim_{index}_{name} :",
                f"    {expression} := by decide",
                "",
            ]
        )
        rows.append(
            {
                "claim_index": index,
                "claim": name,
                "left_exact": str(left),
                "operator": operator,
                "right_exact": str(right),
                "lean_expression": expression,
                "kernel_checked": True,
            }
        )
    lines.extend(
        [
            f"def reflectedClaimCount : Nat := {len(claims)}",
            "",
            "theorem reflected_claim_count :",
            f"    reflectedClaimCount = {len(claims)} := by decide",
            "",
            "end FINPrograms379To393Arithmetic",
            "",
        ]
    )
    ARITHMETIC_SOURCE.write_text("\n".join(lines), encoding="utf-8")
    completed = subprocess.run(
        [str(LEAN_BINARY), str(ARITHMETIC_SOURCE)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
        timeout=180,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr or completed.stdout)
    return (
        {
            "status": (
                "[Proven] exact terminal certificate predicates reflected "
                "through the Lean kernel; [Open] end-to-end algorithm reflection"
            ),
            "kernel_checked_predicate_count": len(claims),
            "krawczyk_containment_predicates": 24,
            "weight_sign_predicates": 7,
            "bernstein_predicates": 2,
            "order_and_gap_predicates": 2,
            "lean_returncode": completed.returncode,
            "arithmetic_source": ARITHMETIC_SOURCE.name,
            "theorem_scope": (
                "Lean checks the exact cross-multiplied integer inequalities "
                "for all P366 Krawczyk containments, all seven signs, the P380 "
                "Bernstein range, weak-duality order, and the sub-1e-7 gap."
            ),
            "boundary": (
                "Python still generates the rational predicates and supplies "
                "the Taylor, fifth-root, interval, and Bernstein algorithms. "
                "This shrinks but does not eliminate the trusted computing base."
            ),
            "new_object": "O121 kernel-reflected terminal certificate",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P380: near-contact exact oscillatory primal-dual enclosure
# ---------------------------------------------------------------------------


def program_380(
    dual_data: dict[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows = []
    for index, (raw, safe, power) in enumerate(
        zip(
            dual_data["raw_bernstein"],
            dual_data["safe_bernstein"],
            dual_data["power"],
        )
    ):
        rows.append(
            {
                "coefficient": index,
                "raw_bernstein_exact": str(raw),
                "safe_bernstein_exact": str(safe),
                "safe_power_exact": str(power),
                "safe_power_decimal": float(power),
            }
        )
    previous = json.loads(
        (ROOT / "FIN_Programs_365_378_Results.json").read_text(encoding="utf-8")
    )["P366"]
    new_enclosure = [
        float(dual_data["dual"].lo),
        float(dual_data["primal"].hi),
    ]
    return (
        {
            "status": "[Proven] exact rational primal-dual gap below 1e-7",
            "bernstein_basis_degree": 11,
            "dyadic_depth": dual_data["depth"],
            "dyadic_cells": 2 ** dual_data["depth"],
            "safe_bernstein_range": [
                float(dual_data["safe_lower"]),
                float(dual_data["safe_upper"]),
            ],
            "certified_enclosure": new_enclosure,
            "certified_gap": float(dual_data["gap"]),
            "gap_below_1e_7": dual_data["gap"] < Fraction(1, 10**7),
            "previous_P366_enclosure": previous["certified_optimum_enclosure"],
            "width_reduction_fraction_from_P366": (
                1.0
                - float(dual_data["gap"])
                / previous["certified_optimum_enclosure_width"]
            ),
            "theorem": (
                "A degree-11 polynomial optimized in the Bernstein basis is "
                "frozen rationally and affine-normalized by its exact depth-14 "
                "dyadic Bernstein range. It lies in [-1,0] on [0,1]. Combined "
                "with the P366 exact seven-atom Krawczyk primal, weak duality "
                "gives the displayed global enclosure."
            ),
            "boundary": (
                "This is a classical signed-moment optimum. Near contact does "
                "not select quantum negativity, information loss, energy, or "
                "a physical coupling."
            ),
            "new_object": "O122 near-contact Bernstein extremal pair",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P381: global integer-distance injectivity
# ---------------------------------------------------------------------------


def strict_kernel(distance: int | float) -> float:
    return math.cos(743.0 * distance / 4000.0 + 13.0 / 80.0) / (
        1.0 + float(distance) ** (9.0 / 5.0)
    )


def program_381() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows = [
        {
            "distance": distance,
            "phase_numerator_over_4000": 743 * distance + 650,
            "kernel_value": strict_kernel(distance),
        }
        for distance in range(65)
    ]
    values = np.array([strict_kernel(distance) for distance in range(4097)])
    rounded_collision_count = len(values) - len(set(values.tolist()))
    return (
        {
            "status": "[Proven] global injectivity on nonnegative integer distances",
            "domain": "d in N_0",
            "finite_diagnostic_range": [0, 4096],
            "finite_double_collision_count": rounded_collision_count,
            "theorem": (
                "Let x_d=(743d+650)/4000 and "
                "a_d=(1+d^(9/5))^(-1). Both are nonzero algebraic numbers "
                "for integer d>=0, and x_d are distinct positive rationals. "
                "If k(d)=k(e), d!=e, then "
                "a_d(e^(ix_d)+e^(-ix_d))-"
                "a_e(e^(ix_e)+e^(-ix_e))=0. The four exponents are distinct "
                "algebraic numbers. Lindemann-Weierstrass linear independence "
                "of their exponentials over the algebraic numbers gives a "
                "contradiction."
            ),
            "consequence": (
                "On every graph with finite integer shortest-path distance, "
                "the P368 strict k-distance quotient is the ordinary distance "
                "set. Strict-kernel-preserving maps are exactly isometric maps."
            ),
            "boundary": (
                "The theorem concerns the frozen strict radial formula and "
                "integer graph distances. It does not identify legacy and "
                "strict kernels or transport physical roles."
            ),
            "new_object": "O123 global strict radial embedding theorem",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P382: exactly one damping-completion naturality square
# ---------------------------------------------------------------------------


def program_382() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows = []
    completed = []
    strict = []
    for distance in range(7):
        legacy_attenuation = 1.0 / (1.0 + 0.01 * distance)
        strict_attenuation = 1.0 / (1.0 + distance ** (9.0 / 5.0))
        completion = (1.0 + 0.01 * distance) / (
            1.0 + distance ** (9.0 / 5.0)
        )
        defect = completion * legacy_attenuation - strict_attenuation
        legacy_phase = math.cos(math.pi * distance / 4.0 + math.pi / 6.0)
        strict_value = strict_kernel(distance)
        completed_value = legacy_phase * strict_attenuation
        rows.append(
            {
                "distance": distance,
                "legacy_attenuation": legacy_attenuation,
                "completion_atom": completion,
                "completed_attenuation": completion * legacy_attenuation,
                "strict_attenuation": strict_attenuation,
                "attenuation_identity_defect": defect,
                "phase_unrepaired_completed_kernel": completed_value,
                "strict_kernel": strict_value,
                "full_kernel_residual": completed_value - strict_value,
            }
        )
        if distance:
            completed.append(completed_value)
            strict.append(strict_value)
    relative_residual = float(
        np.linalg.norm(np.array(completed) - np.array(strict))
        / np.linalg.norm(np.array(strict))
    )
    return (
        {
            "status": (
                "[Proven] one target-defined damping naturality square; "
                "[Refuted] damping atom alone completes the full kernel"
            ),
            "completion_atom": "C_damp(d)=(1+0.01d)/(1+d^(9/5))",
            "exact_identity": (
                "C_damp(d)/(1+0.01d)=1/(1+d^(9/5))"
            ),
            "naturality_theorem": (
                "For an isometric embedding f, pullback of radial matrices "
                "commutes with entrywise multiplication by C_damp because "
                "d_H(fx,fy)=d_G(x,y). This gives one enriched commutative "
                "square for the attenuation layer."
            ),
            "C12_full_kernel_relative_l2_residual_after_damping_only": relative_residual,
            "source_independent": False,
            "boundary": (
                "C_damp reads the strict exponent and scale from its target. "
                "It does not source beta/eta, repair amplitude or phase/frequency, "
                "close the generic bridge, or transfer any legacy physical role."
            ),
            "new_object": "O124 damping-completion natural transformation",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P383: exact coordinate gauge lattice and bounded alias audit
# ---------------------------------------------------------------------------


def program_383(
    strict_a: np.ndarray, rng: np.random.Generator
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    time = float(
        json.loads(
            (ROOT / "FIN_Programs_323_336_Results.json").read_text(
                encoding="utf-8"
            )
        )["P326"]["best_protocols"]["wave"]["nominal_time"]
    )
    target = expm(-1j * time * strict_a)
    rotations, diagonal = p295.givens_decompose_unitary(target)
    count = len(rotations)
    base_parameters = np.zeros(2 * count)
    base = p365.component_transfer(rotations, diagonal, base_parameters)
    rows = []
    period_defects = []
    half_period_defects = []
    loss_defects = []
    for index in range(count):
        period = base_parameters.copy()
        period[count + index] = 2.0 * math.pi
        period_defect = float(
            np.linalg.norm(
                p365.component_transfer(rotations, diagonal, period) - base,
                2,
            )
        )
        half = base_parameters.copy()
        half[count + index] = math.pi
        half_defect = float(
            np.linalg.norm(
                p365.component_transfer(rotations, diagonal, half) - base,
                2,
            )
        )
        loss = base_parameters.copy()
        loss[index] = 0.01
        loss_defect = float(
            np.linalg.norm(
                p365.component_transfer(rotations, diagonal, loss) - base,
                2,
            )
        )
        period_defects.append(period_defect)
        half_period_defects.append(half_defect)
        loss_defects.append(loss_defect)
        rows.append(
            {
                "component": index,
                "two_pi_phase_defect": period_defect,
                "pi_phase_defect": half_defect,
                "loss_shift_0p01_defect": loss_defect,
            }
        )
    combined_defects = []
    for _ in range(32):
        integers = rng.integers(-2, 3, size=count)
        parameters = np.r_[
            np.zeros(count), 2.0 * math.pi * integers.astype(float)
        ]
        defect = float(
            np.linalg.norm(
                p365.component_transfer(rotations, diagonal, parameters) - base,
                2,
            )
        )
        combined_defects.append(defect)
    assert max(period_defects + combined_defects) < 1e-12
    return (
        {
            "status": (
                "[Proven] coordinate phase lattice gauge; "
                "[Strong evidence] no additional tested coordinate aliases; "
                "[Open] global nonlocal alias classification"
            ),
            "component_count": count,
            "parameter_space_before_quotient": "R_+^66 x R^66",
            "certified_coordinate_gauge": "(2*pi*Z)^66 on phase coordinates",
            "canonical_chart": "R_+^66 x (-pi,pi]^66",
            "maximum_single_two_pi_defect": max(period_defects),
            "maximum_combined_lattice_defect": max(combined_defects),
            "minimum_single_pi_defect": min(half_period_defects),
            "minimum_single_loss_shift_defect": min(loss_defects),
            "theorem_scope": (
                "Every component phase enters through complex exponentials, "
                "so independent 2pi shifts are exact gauges. The frozen "
                "coordinatewise pi and positive-loss shifts are not aliases "
                "in the finite audit."
            ),
            "boundary": (
                "Coordinate tests do not exclude compensating, distant, "
                "multi-parameter aliases. Complex transfer tomography and "
                "a global quotient-identifiability theorem remain required."
            ),
            "new_object": "O125 photonic component gauge quotient chart",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P384: physical photonic pilot remains external
# ---------------------------------------------------------------------------


def program_384() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    result, rows = p351.physical_gate(
        "P384", "photonic_hardware_pilot", PHOTONIC_PATH
    )
    result["status"] = (
        "[Blocked by external evidence] no calibrated P384 photonic pilot"
    )
    return result, rows


# ---------------------------------------------------------------------------
# P385: heralded loss/dephasing comb sandwich
# ---------------------------------------------------------------------------


def qubit_phase_state(theta: float, coherence: float, sign: int) -> np.ndarray:
    phase = np.exp(1j * sign * theta)
    return np.array(
        [[0.5, 0.5 * coherence * np.conjugate(phase)],
         [0.5 * coherence * phase, 0.5]],
        dtype=complex,
    )


def half_trace_distance(left: np.ndarray, right: np.ndarray) -> float:
    singular = np.linalg.svd(left - right, compute_uv=False)
    return 0.5 * float(np.sum(singular))


def tensor_power(matrix: np.ndarray, count: int) -> np.ndarray:
    result = np.ones((1, 1), dtype=complex)
    for _ in range(count):
        result = np.kron(result, matrix)
    return result


def product_erasure_lower(
    uses: int, theta: float, survival: float, coherence: float
) -> float:
    plus = qubit_phase_state(theta, coherence, 1)
    minus = qubit_phase_state(theta, coherence, -1)
    total = 0.0
    for successes in range(1, uses + 1):
        weight = (
            math.comb(uses, successes)
            * survival**successes
            * (1.0 - survival) ** (uses - successes)
        )
        distance = half_trace_distance(
            tensor_power(plus, successes),
            tensor_power(minus, successes),
        )
        total += weight * distance
    return total


def program_385() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    previous = json.loads(
        (ROOT / "FIN_Programs_365_378_Results.json").read_text(encoding="utf-8")
    )
    diameter = float(previous["P371"]["relative_generator_spectral_diameter"])
    rows = []
    largest_gap = 0.0
    for uses, survival, coherence, fraction in itertools.product(
        range(1, 5), (0.7, 0.85, 0.95, 1.0), (0.7, 0.9, 1.0), (0.25, 0.5, 0.75, 1.0)
    ):
        threshold = math.pi / (uses * diameter)
        time = fraction * threshold
        theta = time * diameter / 2.0
        one_use = survival * coherence * math.sin(theta)
        adaptive_upper = min(1.0, uses * one_use)
        product_lower = product_erasure_lower(
            uses, theta, survival, coherence
        )
        ghz_lower = (
            survival**uses
            * coherence**uses
            * math.sin(uses * theta)
        )
        lower = max(product_lower, ghz_lower)
        gap = max(0.0, adaptive_upper - lower)
        largest_gap = max(largest_gap, gap)
        if lower > adaptive_upper + 2e-12:
            raise AssertionError((uses, survival, coherence, fraction, lower, adaptive_upper))
        rows.append(
            {
                "uses": uses,
                "survival_probability": survival,
                "coherence_factor": coherence,
                "fraction_of_ideal_threshold": fraction,
                "time_per_use": time,
                "product_erasure_lower": product_lower,
                "ghz_lower": ghz_lower,
                "best_feasible_lower": lower,
                "adaptive_hybrid_upper": adaptive_upper,
                "certificate_gap": gap,
            }
        )
    return (
        {
            "status": (
                "[Proven] analytic upper and feasible lower sandwich for the "
                "declared reduced heralded-loss/dephasing channel; "
                "[Open] exact lossy adaptive optimum"
            ),
            "model": (
                "two extremal relative-generator modes; independent heralded "
                "survival eta; eigenbasis coherence q"
            ),
            "single_use_half_diamond": "eta*q*sin(t*Delta/2)",
            "adaptive_upper": "min(1,n*eta*q*sin(t*Delta/2))",
            "feasible_strategies": [
                "independent extremal probes with erasure-pattern readout",
                "GHZ extremal probe surviving only on the all-click branch",
            ],
            "maximum_grid_certificate_gap": largest_gap,
            "theorem_scope": (
                "The hybrid argument gives the adaptive upper bound. Orthogonal "
                "heralded-erasure blocks give the product lower exactly, and "
                "the all-surviving GHZ block gives the second lower."
            ),
            "boundary": (
                "The bounds do not solve the full twelve-mode lossy comb, "
                "unheralded loss, noncommuting noise, or apparatus calibration."
            ),
            "new_object": "O126 heralded noisy-comb sandwich",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P386: robust clock-tube propagation
# ---------------------------------------------------------------------------


def program_386() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    previous = json.loads(
        (ROOT / "FIN_Programs_365_378_Results.json").read_text(encoding="utf-8")
    )
    diameter = float(previous["P371"]["relative_generator_spectral_diameter"])
    curvature = 0.2
    maximum_gap = 0.2
    epsilon = curvature * maximum_gap**2 / 8.0
    rows = []
    maximum_band = 0.0
    for uses in range(1, 5):
        threshold = math.pi / (uses * diameter)
        for fraction in (0.25, 0.5, 0.75, 1.0):
            nominal = fraction * threshold
            lower_tau = max(0.0, nominal - epsilon)
            upper_tau = min(threshold, nominal + epsilon)
            lower = math.sin(uses * diameter * lower_tau / 2.0)
            upper = math.sin(uses * diameter * upper_tau / 2.0)
            maximum_band = max(maximum_band, upper - lower)
            rows.append(
                {
                    "uses": uses,
                    "fraction_of_first_threshold": fraction,
                    "nominal_tau": nominal,
                    "tau_lower": lower_tau,
                    "tau_upper": upper_tau,
                    "discrimination_lower": lower,
                    "discrimination_upper": upper,
                    "band_width": upper - lower,
                }
            )
    return (
        {
            "status": "[Proven] conditional robust discrimination tube",
            "supplied_curvature_bound": curvature,
            "supplied_maximum_clock_gap": maximum_gap,
            "interpolation_error_bound": epsilon,
            "maximum_reported_discrimination_band_width": maximum_band,
            "theorem": (
                "On the first P371 branch D_n(tau)=sin(n*Delta*tau/2) is "
                "monotone. Therefore a certified tau interval maps exactly "
                "to its endpoint sine interval. The bounded-curvature "
                "interpolation theorem supplies epsilon=M*h^2/8."
            ),
            "boundary": (
                "M, h, the clock anchor, and the relation between laboratory "
                "time and tau are supplied. The theorem propagates calibration "
                "uncertainty; it does not generate a clock or time unit."
            ),
            "new_object": "O127 clock-tube discrimination envelope",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P387: one explicit Jordan-sampling realization functor
# ---------------------------------------------------------------------------


def program_387(
    primal_data: dict[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    centers = primal_data["centers"]
    nodes = np.array(
        [
            float(primal_data["fixed_first"]),
            *[float(value) for value in centers[:5]],
            1.0,
        ]
    )
    weights = np.array([float(value) for value in centers[5:]])
    total_variation = float(np.sum(np.abs(weights)))
    signed_mass = float(np.sum(weights))
    negative_mass = float(-np.sum(weights[weights < 0]))
    probabilities = np.abs(weights) / total_variation
    signs = np.sign(weights)
    moment_targets = np.array(
        [strict_kernel(order) for order in range(12)]
    )
    reconstructed = np.array(
        [
            np.sum(
                probabilities
                * total_variation
                * signs
                * nodes**order
            )
            for order in range(12)
        ]
    )
    residual = float(np.max(np.abs(reconstructed - moment_targets)))
    rows = [
        {
            "atom": index,
            "node": nodes[index],
            "signed_weight": weights[index],
            "sampling_probability": probabilities[index],
            "sign_record": int(signs[index]),
            "importance_score_order_0": total_variation * signs[index],
        }
        for index in range(7)
    ]
    return (
        {
            "status": (
                "[Proven] explicit mathematical operational realization; "
                "[Conditional] physical preparation and calibration"
            ),
            "functor_name": "Jordan Sampling Realization JSR",
            "source_object": "finite signed measure w on [0,1]",
            "preparation": "q_i=|w_i|/||w||_1 with classical sign label",
            "instrument": "record (atom i, node x_i, sign s_i)",
            "score": "||w||_1*s_i*f(x_i)",
            "reconstruction_identity": (
                "E_q[||w||_1*s*f(x)]=sum_i w_i f(x_i)"
            ),
            "free_channel": (
                "positive column-stochastic map T on the signed vector; "
                "||Tw||_1<=||w||_1"
            ),
            "monotone": "N(w)=(||w||_1-sum_i w_i)/2",
            "total_variation": total_variation,
            "signed_mass": signed_mass,
            "negative_mass": negative_mass,
            "negative_label_probability": float(
                np.sum(probabilities[weights < 0])
            ),
            "maximum_center_moment_residual": residual,
            "theorem_scope": (
                "JSR gives preparation, channel class, instrument, score, and "
                "record semantics for the signed-path coordinate. Triangle "
                "inequality proves Jordan negativity is nonincreasing."
            ),
            "boundary": (
                "The sign is an importance-sampling record, not negative "
                "probability, destroyed information, quantum quasiprobability, "
                "or energy. A laboratory realization remains external."
            ),
            "new_object": "O128 Jordan-sampling realization functor",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P388: typed conversion and catalysis matrix
# ---------------------------------------------------------------------------


RESOURCE_NAMES = (
    "signed_path",
    "nontorsion_phase",
    "cross_law",
    "pointing",
    "dimensional_scale",
)


def program_388() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows = []
    cross_refutations = 0
    for source_index, source in enumerate(RESOURCE_NAMES):
        for target_index, target in enumerate(RESOURCE_NAMES):
            countermodel = [0] * len(RESOURCE_NAMES)
            countermodel[source_index] = 1
            if source_index == target_index:
                status = "identity_conversion_accepted"
                counterexample = False
            else:
                status = "universal_positive_cross_conversion_refuted"
                counterexample = True
                cross_refutations += 1
            rows.append(
                {
                    "source": source,
                    "target": target,
                    "status": status,
                    "single_coordinate_countermodel": countermodel,
                    "source_positive": countermodel[source_index] > 0,
                    "target_positive": countermodel[target_index] > 0,
                    "counterexample": counterexample,
                }
            )
    assert cross_refutations == 20
    return (
        {
            "status": (
                "[Proven] no forced positive cross-conversion in the minimal "
                "free resource category; [Open] model-specific conversion laws"
            ),
            "identity_conversions": 5,
            "directed_cross_conversion_implications_refuted": cross_refutations,
            "catalysis": "absent by cancellativity of N^5",
            "reset_maps": (
                "resource-destroying maps to grade zero exist but do not "
                "generate a target resource"
            ),
            "theorem": (
                "For every ordered i!=j, the free generator e_i has positive "
                "i-coordinate and zero j-coordinate, refuting any axiom-free "
                "universal implication i>0 => j>0. If a+c>=b+c in N^5 then "
                "a>=b coordinatewise, so catalysts add no conversions."
            ),
            "boundary": (
                "This is a theorem about the minimal product signature. A "
                "new FIN coupling law or typed operational encoder can add "
                "model-specific conversions and must then be audited separately."
            ),
            "new_object": "O129 resource conversion obstruction matrix",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P389: outer-dilation semigroup and no canonical positive scale
# ---------------------------------------------------------------------------


def program_389() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rho_max = Fraction(12631, 19031)
    rows = []
    maximum_defect = 0.0
    coefficients = np.array([strict_kernel(distance) for distance in range(32)])
    for rho, scale in itertools.product(
        (0.1, 0.25, 0.5, 0.65), (0.25, 0.5, 0.75)
    ):
        if rho >= float(rho_max) or rho * scale >= float(rho_max):
            continue
        direct = coefficients * (rho * scale) ** np.arange(len(coefficients))
        composed = (
            coefficients
            * rho ** np.arange(len(coefficients))
            * scale ** np.arange(len(coefficients))
        )
        defect = float(np.max(np.abs(direct - composed)))
        maximum_defect = max(maximum_defect, defect)
        rows.append(
            {
                "rho": rho,
                "dilation_scale": scale,
                "composed_rho": rho * scale,
                "truncated_semigroup_defect": defect,
                "log_time": -math.log(rho),
                "added_log_time": -math.log(scale),
            }
        )
    return (
        {
            "status": (
                "[Proven] dilation-semigroup law and invariant-selector no-go "
                "for a canonical positive outer scale"
            ),
            "certified_outer_domain_upper": float(rho_max),
            "semigroup": "D_a F_rho = F_(a*rho)",
            "additive_parameter": "t=-log(rho)",
            "maximum_truncated_numerical_defect": maximum_defect,
            "no_go": (
                "A selector invariant under every admissible dilation would "
                "satisfy rho=rho/2, hence rho=0, outside the positive outer "
                "family. The semigroup alone therefore selects no positive rho."
            ),
            "fixed_boundaries": [
                "rho->0 gives the constant coefficient K(0)",
                "rho->1 approaches the undamped generating series boundary",
            ],
            "boundary": (
                "The no-go excludes canonicality from dilation invariance "
                "alone. A new normalization, RG boundary condition, apparatus "
                "resolution, or dimensional axiom could select rho conditionally."
            ),
            "new_object": "O130 outer-scale dilation torsor",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P390--P393: external falsification gates
# ---------------------------------------------------------------------------


def external_gate(
    program: str, kind: str, path: Path, label: str
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    result, rows = p351.physical_gate(program, kind, path)
    result["status"] = f"[Blocked by external evidence] no admitted {label}"
    return result, rows


def run_formal_core() -> dict[str, Any]:
    completed = subprocess.run(
        [str(LEAN_BINARY), str(FORMAL_SOURCE)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr or completed.stdout)
    return {
        "lean_structural_returncode": completed.returncode,
        "lean_arithmetic_returncode": 0,
        "structural_scope": (
            "damping-square naturality, l1 contraction interface, "
            "free-resource cross-conversion obstruction, and scale no-fixed-point"
        ),
        "arithmetic_scope": (
            "terminal exact rational inequalities generated and compiled by P379"
        ),
    }


def summary_rows(results: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "program": f"P{program}",
            "status": results[f"P{program}"]["status"],
            "new_object": results[f"P{program}"].get("new_object", ""),
            "boundary": results[f"P{program}"].get("boundary", ""),
        }
        for program in range(379, 394)
    ]


def make_figures(results: dict[str, Any]) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].bar(
        ["Krawczyk", "Signs", "Bernstein", "Order/gap"],
        [24, 7, 2, 2],
        color=["#176b87", "#64ccc5", "#dafffb", "#d1495b"],
    )
    axes[0].set_ylabel("Lean-reflected predicates")
    axes[0].set_title("P379 arithmetic reflection")
    old_width = (
        results["P380"]["previous_P366_enclosure"][1]
        - results["P380"]["previous_P366_enclosure"][0]
    )
    new_width = results["P380"]["certified_gap"]
    axes[1].bar(
        ["P366", "P380"],
        [old_width, new_width],
        color=["#d1495b", "#176b87"],
    )
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Certified enclosure width")
    axes[1].set_title("P380 contact-gap reduction")
    figure.tight_layout()
    figure.savefig(
        FIGURE_DIR / "p379_p380_reflection_contact.png", dpi=180
    )
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    distances = np.arange(65)
    axes[0].plot(distances, [strict_kernel(int(d)) for d in distances])
    axes[0].axhline(0, color="black", linewidth=0.7)
    axes[0].set_xlabel("integer distance d")
    axes[0].set_ylabel("K_strict(d)")
    axes[0].set_title("P381 globally injective radial samples")
    gauge_rows = results["_P383_rows"]
    axes[1].semilogy(
        [row["component"] for row in gauge_rows],
        [
            max(row["two_pi_phase_defect"], 1e-18)
            for row in gauge_rows
        ],
        label="2pi gauge defect",
    )
    axes[1].semilogy(
        [row["component"] for row in gauge_rows],
        [
            max(row["pi_phase_defect"], 1e-18)
            for row in gauge_rows
        ],
        label="pi shift defect",
    )
    axes[1].set_xlabel("Givens component")
    axes[1].set_title("P383 component phase lattice")
    axes[1].legend()
    figure.tight_layout()
    figure.savefig(
        FIGURE_DIR / "p381_p383_injectivity_gauge.png", dpi=180
    )
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    selected = [
        row
        for row in results["_P385_rows"]
        if row["uses"] == 4
        and row["survival_probability"] == 0.85
        and row["coherence_factor"] == 0.9
    ]
    x = [row["fraction_of_ideal_threshold"] for row in selected]
    axes[0].plot(
        x,
        [row["best_feasible_lower"] for row in selected],
        "o-",
        label="feasible lower",
    )
    axes[0].plot(
        x,
        [row["adaptive_hybrid_upper"] for row in selected],
        "s--",
        label="adaptive upper",
    )
    axes[0].set_xlabel("fraction of ideal threshold")
    axes[0].set_title("P385 noisy-comb sandwich")
    axes[0].legend()
    clock_rows = results["_P386_rows"]
    axes[1].bar(
        [f"n={n}" for n in range(1, 5)],
        [
            max(
                row["band_width"]
                for row in clock_rows
                if row["uses"] == n
            )
            for n in range(1, 5)
        ],
        color="#64ccc5",
    )
    axes[1].set_ylabel("maximum D interval width")
    axes[1].set_title("P386 propagated clock tube")
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p385_p386_loss_clock.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    realization = results["_P387_rows"]
    colors = [
        "#d1495b" if row["sign_record"] < 0 else "#176b87"
        for row in realization
    ]
    axes[0].bar(
        [str(row["atom"]) for row in realization],
        [row["sampling_probability"] for row in realization],
        color=colors,
    )
    axes[0].set_xlabel("signed atom")
    axes[0].set_ylabel("operational probability q")
    axes[0].set_title("P387 Jordan-sampling preparation")
    for rho in (0.1, 0.25, 0.5, 0.65):
        scale = np.linspace(0.1, 1.0, 80)
        axes[1].plot(scale, rho * scale, label=f"rho={rho}")
    axes[1].axhline(
        results["P389"]["certified_outer_domain_upper"],
        color="black",
        linestyle="--",
        label="outer bound",
    )
    axes[1].set_xlabel("dilation a")
    axes[1].set_ylabel("a rho")
    axes[1].set_title("P389 scale orbits")
    axes[1].legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(
        FIGURE_DIR / "p387_p389_realization_scale.png", dpi=180
    )
    plt.close(figure)


def main() -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)
    strict_a, _ = core.strict_operator()
    primal_data = p366_primal_certificate()
    dual_data = p380_dual_certificate(primal_data)

    results: dict[str, Any] = {
        "metadata": {
            "programs": "P379-P393",
            "release": "10.33",
            "seed": SEED,
            "scope": (
                "proof-kernel terminal reflection, near-contact oscillatory "
                "certificate, global distance injectivity, one damping "
                "naturality square, component gauges, noisy combs, robust "
                "clocks, one realization functor, conversion obstructions, "
                "outer-scale no-go, and external falsification gates"
            ),
            "new_theoretical_objects": {
                "O121_kernel_reflected_terminal_certificate": "Lean-checked exact endpoint predicates",
                "O122_near_contact_Bernstein_extremal_pair": "sub-1e-7 signed-moment enclosure",
                "O123_global_strict_radial_embedding_theorem": "integer-distance injectivity",
                "O124_damping_completion_natural_transformation": "one target-defined bridge atom",
                "O125_photonic_component_gauge_quotient_chart": "phase lattice quotient",
                "O126_heralded_noisy_comb_sandwich": "loss/dephasing lower-upper bounds",
                "O127_clock_tube_discrimination_envelope": "robust time-to-discrimination map",
                "O128_Jordan_sampling_realization_functor": "explicit signed-resource operationalization",
                "O129_resource_conversion_obstruction_matrix": "twenty cross-implication countermodels",
                "O130_outer_scale_dilation_torsor": "semigroup canonicality no-go",
            },
        }
    }
    results["P380"], results["_P380_rows"] = program_380(dual_data)
    results["P379"], results["_P379_rows"] = program_379(
        primal_data, dual_data
    )
    results["P381"], results["_P381_rows"] = program_381()
    results["P382"], results["_P382_rows"] = program_382()
    results["P383"], results["_P383_rows"] = program_383(strict_a, rng)
    results["P384"], results["_P384_rows"] = program_384()
    results["P385"], results["_P385_rows"] = program_385()
    results["P386"], results["_P386_rows"] = program_386()
    results["P387"], results["_P387_rows"] = program_387(primal_data)
    results["P388"], results["_P388_rows"] = program_388()
    results["P389"], results["_P389_rows"] = program_389()
    results["P390"], results["_P390_rows"] = external_gate(
        "P390", "qw_independent_holdout", HOLDOUT_PATH, "independent QW hold-out"
    )
    results["P391"], results["_P391_rows"] = external_gate(
        "P391", "dimensional_standards", STANDARDS_PATH, "dimensional standards bundle"
    )
    results["P392"], results["_P392_rows"] = external_gate(
        "P392", "reservoir_process_tomography", RESERVOIR_PATH, "reservoir tomography bundle"
    )
    results["P393"], results["_P393_rows"] = external_gate(
        "P393", "electroweak_blind_test", EW_PATH, "electroweak blind-test bundle"
    )
    results["formal_verification"] = run_formal_core()

    tables = {
        REFLECTION_PATH: results["_P379_rows"],
        CONTACT_PATH: results["_P380_rows"],
        INJECTIVITY_PATH: results["_P381_rows"],
        DAMPING_PATH: results["_P382_rows"],
        GAUGE_PATH: results["_P383_rows"],
        PHOTONIC_PATH: results["_P384_rows"],
        LOSSY_COMB_PATH: results["_P385_rows"],
        CLOCK_PATH: results["_P386_rows"],
        REALIZATION_PATH: results["_P387_rows"],
        CONVERSION_PATH: results["_P388_rows"],
        OUTER_SCALE_PATH: results["_P389_rows"],
        HOLDOUT_PATH: results["_P390_rows"],
        STANDARDS_PATH: results["_P391_rows"],
        RESERVOIR_PATH: results["_P392_rows"],
        EW_PATH: results["_P393_rows"],
    }
    for path, rows in tables.items():
        write_csv(path, rows)
    write_csv(SUMMARY_PATH, summary_rows(results))
    public = {
        key: value for key, value in results.items() if not key.startswith("_")
    }
    RESULTS_PATH.write_text(
        json.dumps(json_ready(public), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    make_figures(results)
    print(RESULTS_PATH)
    for program in range(379, 394):
        print(f"P{program}: {public[f'P{program}']['status']}")


if __name__ == "__main__":
    main()
