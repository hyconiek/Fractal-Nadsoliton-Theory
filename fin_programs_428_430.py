#!/usr/bin/env python3
"""Local analytical/computational batch FIN P428--P430.

P428 repairs and narrows the formal cosine interface.
P429 attempts an exact-rational Krawczyk certificate for the simultaneous
25-equation contact system.
P430 uses the certified contact box, when available, for a contact-aware
global dual feasibility certificate.

No routine in this file consumes external data or emits physical evidence.
"""

from __future__ import annotations

import csv
from fractions import Fraction
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import mpmath as mp
import numpy as np
import matplotlib.pyplot as plt

import fin_programs_351_364 as p351
import fin_programs_394_410 as prior
import fin_programs_411_427 as p411


ROOT = Path(__file__).resolve().parent
PREFIX = "FIN_Programs_428_430"
LEAN = ROOT / ".elan/toolchains/leanprover--lean4---v4.28.0/bin/lean"
LEAN_P428 = ROOT / f"{PREFIX}_P428_Cosine_Reduction.lean"
RESULTS_PATH = ROOT / f"{PREFIX}_Results.json"
SUMMARY_PATH = ROOT / f"{PREFIX}_Summary.csv"
P428_PATH = ROOT / f"{PREFIX}_P428_Cosine.csv"
P429_PATH = ROOT / f"{PREFIX}_P429_Krawczyk.csv"
P430_PATH = ROOT / f"{PREFIX}_P430_Global_Dual.csv"
LEGACY_AUDIT_PATH = ROOT / f"{PREFIX}_LegacyStar_Coupling_Audit.csv"
DIAGRAMS_PATH = ROOT / "DIAGRAMS_KERNEL_TRANSFORMATION.md"
FIGURE_DIR = ROOT / f"{PREFIX}_Figures"
FIGURE_PATH = FIGURE_DIR / "p428_p430_exact_closure_and_legacy_context.png"


def json_ready(value: Any) -> Any:
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, mp.mpf):
        return mp.nstr(value, 50)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("status\nno_rows\n", encoding="utf-8")
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
            writer.writerow({
                key: json.dumps(json_ready(row.get(key)), ensure_ascii=False)
                if isinstance(row.get(key), (dict, list, tuple))
                else row.get(key, "")
                for key in keys
            })


def mp_decimal(value: mp.mpf, digits: int = 100) -> str:
    return mp.nstr(value, n=digits, min_fixed=-100000, max_fixed=100000)


def mp_to_fraction(value: mp.mpf, digits: int = 90) -> Fraction:
    return Fraction(mp_decimal(value, digits))


def ri_abs_upper(value: p351.RI) -> Fraction:
    return max(abs(value.lo), abs(value.hi))


def ri_hull(values: list[p351.RI]) -> p351.RI:
    return p351.RI(min(value.lo for value in values), max(value.hi for value in values))


def legacy_star_coupling_audit() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Audit the historical coupling diagram used to reconstruct Legacy*.

    This is a source/provenance and algebra check, not a new physics program.
    It keeps the reconstructed historical kernel visible in a strict-kernel
    batch without silently transferring its physical-role prose.
    """

    text_source = DIAGRAMS_PATH.read_text(encoding="utf-8")
    source_checks = {
        "four_factor_product_present": (
            "K_total = K_geo × K_res × (1+0.2K_tors) × K_topo" in text_source
        ),
        "effective_legacy_formula_present": (
            "K(d) = α·cos(ωd+φ)/(1+β·d)" in text_source
        ),
        "resonance_feedback_present": (
            "α_res_eff = α_res_base × exp(-0.5·α_geo)" in text_source
        ),
        "topology_torsion_feedback_present": (
            "β_topo_eff = β_topo_base × (1 + 0.5|K_tors|)" in text_source
        ),
    }
    if not all(source_checks.values()):
        raise RuntimeError(f"Legacy coupling source check failed: {source_checks}")

    phase_values = {
        distance: math.cos(math.pi * distance / 4 + math.pi / 6)
        for distance in (2, 5, 8, 11)
    }
    rows = [
        {
            "item": "K_geo",
            "declared_role": "viscosity / damping",
            "diagram_form": "exp(-alpha*d)",
            "audited_status": "historical mechanism label; not a physical source theorem",
        },
        {
            "item": "K_res",
            "declared_role": "resonance / phase synchronization",
            "diagram_form": "1+alpha_res*similarity",
            "audited_status": "historical mechanism label; receiver-like unless a state law is supplied",
        },
        {
            "item": "K_tors",
            "declared_role": "oscillatory currents",
            "diagram_form": "cos(pi*d/4+pi/6)",
            "audited_status": "exact scalar carrier; no orientation selector is exported",
        },
        {
            "item": "K_topo",
            "declared_role": "topological/path coupling",
            "diagram_form": "path sum -> 1/(1+beta*d)",
            "audited_status": "conditional reconstruction; the displayed exponent derivation is algebraically wrong",
        },
        {
            "item": "interdependence",
            "declared_role": "cross-modulation of effective resonance and topology parameters",
            "diagram_form": "alpha_res_eff(alpha_geo), beta_topo_eff(K_tors)",
            "audited_status": "self-coupling intuition only: no typed self-map, feedback state equation, or fixed-point law",
        },
        {
            "item": "LegacyStar",
            "declared_role": "reconstructed historical effective kernel",
            "diagram_form": "A*cos(pi*d/4+pi/6)/(1+beta*d)",
            "audited_status": "accepted kernel class; A and beta remain free/frozen; intermediate relative to strict",
        },
    ]
    return ({
        "status": (
            "[Proven audit] the diagram contains four coupling mechanisms and two "
            "cross-modulation formulas; [Refuted] its displayed path exponent and "
            "integer-node derivations; [Conditional] these formulas are a self-"
            "coupling intuition, not a typed dynamical self-coupling theorem"
        ),
        "source_checks": source_checks,
        "legacy_star_definition": "A*cos(pi*d/4+pi/6)/(1+beta*d), A>0, beta>0",
        "relation_to_binding_legacy": (
            "the fixed-phase subclass of K_legacy_ont with alpha_geo=A, "
            "omega=pi/4, phi=pi/6, beta_tors=beta"
        ),
        "path_sum_exponent_claim": 1.6 + (-0.6),
        "path_sum_exponent_required_for_inverse_linear_tail": -2.6,
        "alleged_integer_node_cosine_values": phase_values,
        "exact_real_zero_set": "d=4/3+4n",
        "self_coupling_type_boundary": (
            "No map S:X->X, no state variable, no evolution equation, no feedback "
            "closure, and no fixed-point condition is defined in the diagram."
        ),
        "physical_evidence_boundary": (
            "Historical Wilson-loop and accuracy statements in the local diagram "
            "are repository claims, not admitted external laboratory evidence."
        ),
        "selector_discharged": False,
        "physical_role_transfer_exported": False,
    }, rows)


def make_checkpoint_figure(
    rows_428: list[dict[str, Any]], payload: dict[str, Any]
) -> None:
    """Create a compact visual audit of the exact batch and legacy context."""

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    centers: list[Fraction] = payload["centers"]
    nodes = np.array([float(value) for value in centers[:6]] + [1.0])
    weights = np.array([float(value) for value in centers[6:13]])
    coefficients = np.array([float(value) for value in centers[13:25]])
    grid = np.linspace(0.0, 1.0, 4001)
    polynomial = np.polynomial.polynomial.polyval(grid, coefficients)

    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.6), constrained_layout=True)

    ax = axes[0, 0]
    ax.semilogy(
        [row["distance"] for row in rows_428],
        [row["width"] for row in rows_428],
        "o-",
        color="#1f5a99",
    )
    ax.axhline(1e-30, color="#9b2c2c", linestyle="--", linewidth=1.1)
    ax.set_title("P428: rational cosine enclosure widths")
    ax.set_xlabel("distance d")
    ax.set_ylabel("upper - lower")
    ax.grid(alpha=0.25)

    ax = axes[0, 1]
    colors = ["#9b2c2c" if value < 0 else "#247a45" for value in weights]
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.scatter(nodes, weights, c=colors, s=55, zorder=3)
    for x_value, weight in zip(nodes, weights):
        ax.vlines(x_value, 0, weight, color="#555555", linewidth=1.0)
    ax.set_title("P429: certified seven-atom sign pattern")
    ax.set_xlabel("contact node")
    ax.set_ylabel("signed weight")
    ax.grid(alpha=0.2)

    ax = axes[1, 0]
    ax.fill_between(grid, -1.0, 0.0, color="#dcebdc", alpha=0.75, label="feasible strip")
    ax.plot(grid, polynomial, color="#3c277a", linewidth=1.6, label="dual p(x), midpoint")
    ax.scatter(nodes, np.polynomial.polynomial.polyval(nodes, coefficients), color="#c05a18", s=24, zorder=4)
    ax.set_ylim(-1.08, 0.08)
    ax.set_title("P430: globally certified dual feasibility")
    ax.set_xlabel("x")
    ax.set_ylabel("p(x)")
    ax.legend(loc="lower center", fontsize=8)
    ax.grid(alpha=0.2)

    ax = axes[1, 1]
    ax.axis("off")
    factors = [
        (0.08, 0.72, "$K_{geo}$\ndamping"),
        (0.31, 0.72, "$K_{res}$\nresonance"),
        (0.54, 0.72, "$K_{tors}$\noscillation"),
        (0.77, 0.72, "$K_{topo}$\npath sum"),
    ]
    for x_value, y_value, label in factors:
        ax.text(
            x_value,
            y_value,
            label,
            transform=ax.transAxes,
            ha="center",
            va="center",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#edf2f7", edgecolor="#56708a"),
            fontsize=9,
        )
        ax.annotate(
            "",
            xy=(0.5, 0.35),
            xytext=(x_value, y_value - 0.08),
            xycoords=ax.transAxes,
            textcoords=ax.transAxes,
            arrowprops=dict(arrowstyle="->", color="#56708a", linewidth=1.1),
        )
    ax.text(
        0.5,
        0.28,
        r"$K^*_{legacy}(d)=A\cos(\pi d/4+\pi/6)/(1+\beta d)$",
        transform=ax.transAxes,
        ha="center",
        va="center",
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#fff4dc", edgecolor="#b07b1d"),
        fontsize=10,
    )
    ax.text(
        0.5,
        0.08,
        "Reconstructed class, not a typed self-coupling law;\nA and beta remain unsourced parameters.",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=8.5,
        color="#7a2b2b",
    )
    ax.set_title("Legacy* coupling provenance and boundary")

    fig.suptitle("FIN local checkpoint P428–P430: exact closure and guarded legacy context", fontsize=14)
    fig.savefig(FIGURE_PATH, dpi=210)
    plt.close(fig)


# ---------------------------------------------------------------------------
# P428: type-correct formal dependency reduction
# ---------------------------------------------------------------------------


def program_428() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    completed = subprocess.run(
        [str(LEAN), LEAN_P428.name],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr or completed.stdout)

    library_root = ROOT / ".elan/toolchains/leanprover--lean4---v4.28.0/lib/lean"
    real_candidates = [
        path for path in library_root.rglob("*")
        if path.is_file()
        and path.suffix in {".lean", ".olean"}
        and any(token in path.name.lower() for token in ("real", "trig", "cosine"))
    ]
    mathlib_present = (ROOT / "Mathlib.lean").exists() or (ROOT / "Mathlib").is_dir()

    rows: list[dict[str, Any]] = []
    maximum_width = Fraction(0)
    for distance in range(12):
        angle = Fraction(743 * distance + 650, 4000)
        lower = p411.rational_cos_partial(angle, 21)
        upper = p411.rational_cos_partial(angle, 20)
        maximum_width = max(maximum_width, upper - lower)
        rows.append({
            "distance": distance,
            "angle_exact": str(angle),
            "nonnegative": angle >= 0,
            "square_below_12": angle * angle < 12,
            "lower_exact": str(lower),
            "upper_exact": str(upper),
            "width": float(upper - lower),
        })

    return ({
        "status": (
            "[Proven] FIN-specific rational side conditions and type-correct "
            "cut-interface reduction compile in Lean; [Refuted] Rat->Rat is a "
            "type-adequate interface for standard real cosine; [Blocked locally] "
            "the standard real-analysis cosine theorem is absent"
        ),
        "lean_returncode": completed.returncode,
        "angle_count": 12,
        "maximum_rational_width": float(maximum_width),
        "mathlib_present": mathlib_present,
        "local_real_trig_candidate_files": len(real_candidates),
        "old_interface": "Rat -> Rat",
        "old_interface_type_adequate_for_real_cosine": False,
        "replacement_interface": "Rat -> CutValue(lowerBound, upperBound, compatibility)",
        "formal_reduction": (
            "All twelve FIN bounds follow from one generic standard alternating-"
            "cosine provider on 0<=x and x^2<12; every FIN-specific side "
            "condition is proof-kernel checked."
        ),
        "boundary": (
            "The local Lean/Std toolchain contains no standard real-number "
            "trigonometric analysis library. The theorem relating the power "
            "series to the usual real cosine remains an external mathematical "
            "provider, not a FIN axiom and not an internally compiled theorem."
        ),
        "new_object": "O155 Rational-Cut Cosine Provider Interface",
    }, rows)


# ---------------------------------------------------------------------------
# P429: high-precision root plus exact-rational Krawczyk inclusion
# ---------------------------------------------------------------------------


def strict_moments_mp() -> list[mp.mpf]:
    return [
        mp.cos(mp.mpf(743) * order / 4000 + mp.mpf(13) / 80)
        / (1 + mp.power(order, mp.mpf(9) / 5))
        for order in range(12)
    ]


def kkt_system_mp(vector: list[mp.mpf], moments: list[mp.mpf]) -> list[mp.mpf]:
    nodes = vector[:6] + [mp.mpf(1)]
    weights = vector[6:13]
    coefficients = vector[13:25]
    targets = [mp.mpf(-1), 0, 0, 0, 0, mp.mpf(-1), 0]
    residuals: list[mp.mpf] = []
    for order in range(12):
        residuals.append(mp.fsum(
            weights[index] * nodes[index] ** order for index in range(7)
        ) - moments[order])
    for index, node in enumerate(nodes):
        residuals.append(mp.fsum(
            coefficients[power] * node**power for power in range(12)
        ) - targets[index])
    for index in range(6):
        node = nodes[index]
        residuals.append(mp.fsum(
            power * coefficients[power] * node ** (power - 1)
            for power in range(1, 12)
        ))
    return residuals


def kkt_jacobian_mp(vector: list[mp.mpf]) -> mp.matrix:
    nodes = vector[:6] + [mp.mpf(1)]
    weights = vector[6:13]
    coefficients = vector[13:25]
    jacobian = mp.matrix(25, 25)
    for order in range(12):
        for index in range(6):
            jacobian[order, index] = (
                weights[index] * order * nodes[index] ** (order - 1)
                if order else 0
            )
        for index in range(7):
            jacobian[order, 6 + index] = nodes[index] ** order
    for index, node in enumerate(nodes):
        row = 12 + index
        if index < 6:
            jacobian[row, index] = mp.fsum(
                power * coefficients[power] * node ** (power - 1)
                for power in range(1, 12)
            )
        for power in range(12):
            jacobian[row, 13 + power] = node**power
    for index in range(6):
        row = 19 + index
        node = nodes[index]
        jacobian[row, index] = mp.fsum(
            power * (power - 1) * coefficients[power] * node ** (power - 2)
            for power in range(2, 12)
        )
        for power in range(1, 12):
            jacobian[row, 13 + power] = power * node ** (power - 1)
    return jacobian


def refine_kkt_root() -> tuple[list[mp.mpf], dict[str, Any]]:
    mp.mp.dps = 180
    vector = [
        *map(mp.mpf, prior.KKT_NODES[:6]),
        *map(mp.mpf, prior.KKT_WEIGHTS),
        *map(mp.mpf, prior.KKT_POWER),
    ]
    moments = strict_moments_mp()
    history: list[str] = []
    for _ in range(20):
        residual = kkt_system_mp(vector, moments)
        norm = max(abs(value) for value in residual)
        history.append(mp.nstr(norm, 12))
        if norm < mp.mpf("1e-145"):
            break
        jacobian = kkt_jacobian_mp(vector)
        delta = mp.lu_solve(jacobian, mp.matrix([-value for value in residual]))
        vector = [vector[index] + delta[index] for index in range(25)]
    residual = kkt_system_mp(vector, moments)
    jacobian_np = np.array(
        [[float(value) for value in row] for row in kkt_jacobian_mp(vector).tolist()],
        dtype=float,
    )
    singular = np.linalg.svd(jacobian_np, compute_uv=False)
    return vector, {
        "newton_iterations": len(history),
        "residual_history": history,
        "maximum_residual": mp.nstr(max(abs(value) for value in residual), 30),
        "smallest_double_singular_value": float(singular[-1]),
        "double_condition_number": float(singular[0] / singular[-1]),
    }


def kkt_system_ri(variables: list[p351.RI], moments: list[p351.RI]) -> list[p351.RI]:
    nodes = variables[:6] + [p351.RI.point(1)]
    weights = variables[6:13]
    coefficients = variables[13:25]
    targets = [Fraction(-1), 0, 0, 0, 0, Fraction(-1), 0]
    result: list[p351.RI] = []
    for order in range(12):
        result.append(p351.interval_sum([
            weights[index] * nodes[index].power(order) for index in range(7)
        ]) - moments[order])
    for index, node in enumerate(nodes):
        result.append(p351.interval_sum([
            coefficients[power] * node.power(power) for power in range(12)
        ]) - p351.RI.point(targets[index]))
    for index in range(6):
        node = nodes[index]
        result.append(p351.interval_sum([
            coefficients[power]
            * node.power(power - 1)
            * p351.RI.point(power)
            for power in range(1, 12)
        ]))
    return result


def kkt_jacobian_ri(variables: list[p351.RI]) -> list[list[p351.RI]]:
    zero = p351.RI.point(0)
    nodes = variables[:6] + [p351.RI.point(1)]
    weights = variables[6:13]
    coefficients = variables[13:25]
    rows = [[zero for _ in range(25)] for _ in range(25)]
    for order in range(12):
        for index in range(6):
            rows[order][index] = (
                weights[index]
                * nodes[index].power(order - 1)
                * p351.RI.point(order)
                if order else zero
            )
        for index in range(7):
            rows[order][6 + index] = nodes[index].power(order)
    for index, node in enumerate(nodes):
        row = 12 + index
        if index < 6:
            rows[row][index] = p351.interval_sum([
                coefficients[power]
                * node.power(power - 1)
                * p351.RI.point(power)
                for power in range(1, 12)
            ])
        for power in range(12):
            rows[row][13 + power] = node.power(power)
    for index in range(6):
        row = 19 + index
        node = nodes[index]
        rows[row][index] = p351.interval_sum([
            coefficients[power]
            * node.power(power - 2)
            * p351.RI.point(power * (power - 1))
            for power in range(2, 12)
        ])
        for power in range(1, 12):
            rows[row][13 + power] = (
                node.power(power - 1) * p351.RI.point(power)
            )
    return rows


def exact_krawczyk(
    refined: list[mp.mpf],
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any] | None]:
    centers = [mp_to_fraction(value, 100) for value in refined]
    point = [p351.RI.point(value) for value in centers]
    moments = [p351.oscillatory_moment_interval(order) for order in range(12)]
    # A Krawczyk preconditioner need not be the exact inverse of the rational
    # midpoint Jacobian: any fixed matrix C is valid.  Computing that inverse
    # by Fraction Gaussian elimination causes avoidable expression swell for
    # this ill-conditioned 25-dimensional system.  We therefore invert the
    # high-precision midpoint Jacobian numerically once, round every entry to
    # an explicitly stored rational decimal, and perform every subsequent
    # inclusion operation exactly over Fraction intervals.  Rounding cannot
    # invalidate the proof because the rounded C itself is the matrix used in
    # the exact Krawczyk map.
    inverse_mp = kkt_jacobian_mp(refined) ** -1
    inverse = [
        [mp_to_fraction(inverse_mp[row, column], 70) for column in range(25)]
        for row in range(25)
    ]
    point_system = kkt_system_ri(point, moments)
    correction = p351.matrix_vector_interval(inverse, point_system)

    attempts: list[dict[str, Any]] = []
    success_payload: dict[str, Any] | None = None
    for exponent in (20, 24, 28, 30, 32, 34, 36, 38, 40, 42):
        relative = Fraction(1, 10**exponent)
        radii = [relative * max(Fraction(1), abs(center)) for center in centers]
        box = [
            p351.RI(center - radius, center + radius)
            for center, radius in zip(centers, radii)
        ]
        jacobian_box = kkt_jacobian_ri(box)
        cj = p351.interval_matrix_product(inverse, jacobian_box)
        delta = [p351.RI(-radius, radius) for radius in radii]
        contraction_rows = [
            sum(
                ri_abs_upper(
                    p351.RI.point(int(row == column)) - cj[row][column]
                )
                * radii[column]
                / radii[row]
                for column in range(25)
            )
            for row in range(25)
        ]
        contraction_upper = max(contraction_rows)
        base = [point[index] - correction[index] for index in range(25)]
        image: list[p351.RI] = []
        for row in range(25):
            remainder = p351.interval_sum([
                (p351.RI.point(int(row == column)) - cj[row][column])
                * delta[column]
                for column in range(25)
            ])
            image.append(base[row] + remainder)
        inside = [image[index].strictly_inside(box[index]) for index in range(25)]
        width_ratios = [
            float(image[index].width / box[index].width)
            for index in range(25)
        ]
        margins = [
            min(
                image[index].lo - box[index].lo,
                box[index].hi - image[index].hi,
            ) / radii[index]
            for index in range(25)
        ]
        attempt = {
            "relative_radius_exponent": exponent,
            "inside_variables": sum(inside),
            "all_inside": all(inside),
            "maximum_image_to_box_width_ratio": max(width_ratios),
            "minimum_normalized_inclusion_margin": float(min(margins)),
            "weighted_infinity_contraction_upper": float(contraction_upper),
            "strict_contraction": contraction_upper < 1,
        }
        attempts.append(attempt)
        if all(inside) and contraction_upper < 1:
            signs = [
                1 if cell.lo > 0 else -1 if cell.hi < 0 else 0
                for cell in image[6:13]
            ]
            ordered = all(
                image[index].hi < image[index + 1].lo for index in range(5)
            ) and image[5].hi < 1
            success_payload = {
                "centers": centers,
                "radii": radii,
                "box": box,
                "image": image,
                "inverse": inverse,
                "moments": moments,
                "weight_signs": signs,
                "node_order_certified": ordered,
                "relative_radius_exponent": exponent,
                "weighted_infinity_contraction_upper": contraction_upper,
            }
            break

    rows: list[dict[str, Any]] = []
    for attempt in attempts:
        rows.append({"row_type": "radius_attempt", **attempt})
    if success_payload is not None:
        labels = (
            [f"x{index}" for index in range(6)]
            + [f"w{index}" for index in range(7)]
            + [f"c{index}" for index in range(12)]
        )
        for label, cell, image in zip(
            labels, success_payload["box"], success_payload["image"]
        ):
            rows.append({
                "row_type": "certified_variable",
                "variable": label,
                "box_lower": str(cell.lo),
                "box_upper": str(cell.hi),
                "image_lower": str(image.lo),
                "image_upper": str(image.hi),
                "strictly_inside": image.strictly_inside(cell),
                "box_width_float": float(cell.width),
                "image_width_float": float(image.width),
            })

    result = {
        "status": (
            "[Computer-assisted proof] unique seven-contact KKT zero in the "
            "declared exact-rational box"
            if success_payload is not None
            else "[Open] exact-rational Krawczyk inclusion not obtained on the declared radius ladder"
        ),
        "attempts": len(attempts),
        "successful_relative_radius_exponent": (
            success_payload["relative_radius_exponent"]
            if success_payload is not None else None
        ),
        "all_25_variables_strictly_included": success_payload is not None,
        "weight_signs": (
            success_payload["weight_signs"] if success_payload is not None else None
        ),
        "node_order_certified": (
            success_payload["node_order_certified"] if success_payload is not None else False
        ),
        "weighted_infinity_contraction_upper": (
            float(success_payload["weighted_infinity_contraction_upper"])
            if success_payload is not None else None
        ),
        "moment_enclosure": (
            "exact Fraction arithmetic; rational Taylor/Lagrange cosine bounds "
            "and exact fifth-root brackets"
        ),
        "preconditioner": (
            "70-decimal rational rounding of a 180-digit numerical inverse; "
            "the rounded matrix is used exactly in the Krawczyk map"
        ),
        "boundary": (
            "The certificate concerns the 25-equation contact system and this "
            "isolating box. Global dual feasibility and global optimality are "
            "separate P430 obligations."
        ),
        "new_object": "O156 Exact-Rational Seven-Contact KKT Isolating Box",
    }
    return result, rows, success_payload


def program_429() -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any] | None]:
    refined, refinement = refine_kkt_root()
    result, rows, payload = exact_krawczyk(refined)
    objective_box = (
        -payload["box"][6] - payload["box"][11]
        if payload is not None else None
    )
    result.update({
        "refinement": refinement,
        "objective_midpoint": float(-refined[6] - refined[11]),
        "objective_interval_lower": (
            float(objective_box.lo) if objective_box is not None else None
        ),
        "objective_interval_upper": (
            float(objective_box.hi) if objective_box is not None else None
        ),
        "objective_interval_width": (
            float(objective_box.width) if objective_box is not None else None
        ),
    })
    return result, rows, payload


# ---------------------------------------------------------------------------
# P430: contact-aware global polynomial feasibility
# ---------------------------------------------------------------------------


def derivative_coefficients(
    coefficients: list[p351.RI], order: int = 1
) -> list[p351.RI]:
    result = coefficients
    for _ in range(order):
        result = [
            result[index] * p351.RI.point(index)
            for index in range(1, len(result))
        ]
    return result


def affine_power_interval(
    coefficients: list[p351.RI], left: Fraction, right: Fraction
) -> list[p351.RI]:
    degree = len(coefficients) - 1
    width = right - left
    result = [p351.RI.point(0) for _ in range(degree + 1)]
    for power, coefficient in enumerate(coefficients):
        for t_power in range(power + 1):
            factor = (
                Fraction(math.comb(power, t_power))
                * left ** (power - t_power)
                * width**t_power
            )
            result[t_power] = result[t_power] + coefficient.scale(factor)
    return result


def power_to_bernstein_interval(power: list[p351.RI]) -> list[p351.RI]:
    degree = len(power) - 1
    result: list[p351.RI] = []
    for index in range(degree + 1):
        result.append(p351.interval_sum([
            power[power_index].scale(
                Fraction(math.comb(index, power_index), math.comb(degree, power_index))
            )
            for power_index in range(index + 1)
        ]))
    return result


def bernstein_range_interval(
    coefficients: list[p351.RI], left: Fraction, right: Fraction
) -> p351.RI:
    local_power = affine_power_interval(coefficients, left, right)
    bernstein = power_to_bernstein_interval(local_power)
    return ri_hull(bernstein)


def adaptive_gap_certificate(
    coefficients: list[p351.RI],
    left: Fraction,
    right: Fraction,
    depth: int,
    maximum_depth: int,
    rows: list[dict[str, Any]],
) -> bool:
    bound = bernstein_range_interval(coefficients, left, right)
    safe = bound.lo > -1 and bound.hi < 0
    if safe:
        rows.append({
            "row_type": "noncontact_bernstein_cell",
            "left": float(left),
            "right": float(right),
            "lower": float(bound.lo),
            "upper": float(bound.hi),
            "depth": depth,
            "safe": True,
        })
        return True
    if depth >= maximum_depth:
        rows.append({
            "row_type": "noncontact_bernstein_cell",
            "left": float(left),
            "right": float(right),
            "lower": float(bound.lo),
            "upper": float(bound.hi),
            "depth": depth,
            "safe": False,
        })
        return False
    middle = (left + right) / 2
    return (
        adaptive_gap_certificate(
            coefficients, left, middle, depth + 1, maximum_depth, rows
        )
        and adaptive_gap_certificate(
            coefficients, middle, right, depth + 1, maximum_depth, rows
        )
    )


def program_430(
    payload: dict[str, Any] | None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if payload is None:
        return ({
            "status": "[Blocked by P429] no certified KKT box available for contact-aware global proof",
            "global_feasibility_certified": False,
            "boundary": "A floating-point polynomial is not substituted for the missing interval root object.",
            "new_object": None,
        }, [{"row_type": "dependency", "P429_certificate_available": False}])

    # Use the simple decimal isolating box, not the much more complicated
    # exact endpoints of its Krawczyk image.  The inclusion K(X) subset int(X)
    # already proves that the unique root lies in X; retaining X here avoids
    # catastrophic rational-expression growth without weakening the bound.
    root_box: list[p351.RI] = payload["box"]
    node_boxes = root_box[:6] + [p351.RI.point(1)]
    coefficients = root_box[13:25]
    second = derivative_coefficients(coefficients, 2)
    first = derivative_coefficients(coefficients, 1)
    neighborhood_halfwidth = Fraction(1, 10**4)
    rows: list[dict[str, Any]] = []
    contact_safe: list[bool] = []

    # Six stationary contacts: negative minima at 0 and 5, zero maxima at 1..4.
    neighborhoods: list[tuple[Fraction, Fraction]] = []
    for index, node in enumerate(node_boxes[:6]):
        left = max(Fraction(0), node.lo - neighborhood_halfwidth)
        right = min(Fraction(1), node.hi + neighborhood_halfwidth)
        neighborhoods.append((left, right))
        second_bound = bernstein_range_interval(second, left, right)
        value_bound = bernstein_range_interval(coefficients, left, right)
        if index in (0, 5):
            curvature = second_bound.lo > 0
            opposite = value_bound.hi < 0
            kind = "negative_contact_convex_minimum"
        else:
            curvature = second_bound.hi < 0
            opposite = value_bound.lo > -1
            kind = "zero_contact_concave_maximum"
        safe = curvature and opposite
        contact_safe.append(safe)
        rows.append({
            "row_type": kind,
            "contact": index,
            "left": float(left),
            "right": float(right),
            "second_derivative_lower": float(second_bound.lo),
            "second_derivative_upper": float(second_bound.hi),
            "direct_value_lower": float(value_bound.lo),
            "direct_value_upper": float(value_bound.hi),
            "safe": safe,
        })

    # Endpoint x=1 is a zero contact.  Positive derivative to its left proves p<=0.
    endpoint_left = Fraction(1) - neighborhood_halfwidth
    endpoint_right = Fraction(1)
    endpoint_derivative = bernstein_range_interval(first, endpoint_left, endpoint_right)
    endpoint_value = bernstein_range_interval(coefficients, endpoint_left, endpoint_right)
    endpoint_safe = endpoint_derivative.lo > 0 and endpoint_value.lo > -1
    rows.append({
        "row_type": "endpoint_zero_monotone",
        "contact": 6,
        "left": float(endpoint_left),
        "right": 1.0,
        "derivative_lower": float(endpoint_derivative.lo),
        "derivative_upper": float(endpoint_derivative.hi),
        "direct_value_lower": float(endpoint_value.lo),
        "direct_value_upper": float(endpoint_value.hi),
        "safe": endpoint_safe,
    })

    excluded = neighborhoods + [(endpoint_left, endpoint_right)]
    excluded.sort()
    gaps: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for left, right in excluded:
        if cursor < left:
            gaps.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < 1:
        gaps.append((cursor, Fraction(1)))

    gap_safe = []
    for left, right in gaps:
        gap_safe.append(adaptive_gap_certificate(
            coefficients, left, right, 0, 18, rows
        ))

    global_safe = all(contact_safe) and endpoint_safe and all(gap_safe)
    return ({
        "status": (
            "[Computer-assisted proof] the unique P429 dual obeys -1<=p<=0 on [0,1]"
            if global_safe
            else "[Open] contact-aware interval decomposition did not certify every global cell"
        ),
        "global_feasibility_certified": global_safe,
        "stationary_contacts_certified": sum(contact_safe),
        "stationary_contacts_total": 6,
        "endpoint_certified": endpoint_safe,
        "noncontact_gaps_certified": sum(gap_safe),
        "noncontact_gaps_total": len(gap_safe),
        "bernstein_terminal_cells": sum(
            row["row_type"] == "noncontact_bernstein_cell" for row in rows
        ),
        "maximum_subdivision_depth": max(
            (row.get("depth", 0) for row in rows), default=0
        ),
        "theorem_scope": (
            "Krawczyk uniqueness supplies exact contact equations; interval "
            "curvature handles the six stationary contacts, interval monotonicity "
            "handles x=1, and Bernstein ranges handle the complement."
        ),
        "boundary": (
            "This certifies dual feasibility for the P429 contact solution. "
            "Primal moment feasibility is part of the same KKT zero, but global "
            "optimality also uses the signed-measure duality theorem inherited "
            "from P337--P380. It creates no physical interpretation."
        ),
        "new_object": "O157 Contact-Aware Global Dual Certificate" if global_safe else None,
    }, rows)


def main() -> None:
    results: dict[str, Any] = {
        "metadata": {
            "programs": "P428-P430",
            "checkpoint": "P428-P430",
            "date": "2026-08-01",
            "local_only": True,
            "external_physical_evidence": False,
            "kernel_split_preserved": True,
            "selector_discharged": False,
            "dimensional_source_exported": False,
        }
    }
    results["legacy_star_coupling_context"], legacy_rows = legacy_star_coupling_audit()
    results["P428"], rows_428 = program_428()
    results["P429"], rows_429, payload = program_429()
    results["P430"], rows_430 = program_430(payload)
    if payload is not None:
        make_checkpoint_figure(rows_428, payload)
    write_csv(P428_PATH, rows_428)
    write_csv(P429_PATH, rows_429)
    write_csv(P430_PATH, rows_430)
    write_csv(LEGACY_AUDIT_PATH, legacy_rows)
    write_csv(SUMMARY_PATH, [
        {"program": name, "status": results[name]["status"]}
        for name in ("P428", "P429", "P430")
    ])
    RESULTS_PATH.write_text(
        json.dumps(json_ready(results), indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        name: results[name]["status"] for name in ("P428", "P429", "P430")
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
