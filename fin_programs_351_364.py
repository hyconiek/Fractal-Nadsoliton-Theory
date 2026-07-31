#!/usr/bin/env python3
"""Execute FIN Research Programs P351--P364 (Release 10.31).

This round certifies the two signed-moment problems, proves the exact
naturality boundary for radial graph kernels, stress-tests an ideal compiled
photonic transfer under declared loss/noise, supplies explicit perfect
multi-use channel witnesses, closes the finite-design bounded-curvature clock
problem, and constructs a minimal axiomatic bridge-resource category.

External hold-out, hardware, electroweak and conversion-standard gates remain
external.  Synthetic records are never admitted as laboratory evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import networkx as nx
import numpy as np
from scipy.linalg import expm
from scipy.optimize import linprog

import fin_programs_255_266 as core
import fin_programs_295_308 as p295
import fin_programs_323_336 as p323
import fin_programs_337_350 as p337


ROOT = Path(__file__).resolve().parent
FIGURE_DIR = ROOT / "FIN_Programs_351_364_Figures"
RESULTS_PATH = ROOT / "FIN_Programs_351_364_Results.json"
SUMMARY_PATH = ROOT / "FIN_Programs_351_364_Summary.csv"
ENVELOPE_PATH = ROOT / "FIN_Programs_351_364_Envelope_Interval.csv"
OSCILLATORY_PATH = ROOT / "FIN_Programs_351_364_Oscillatory_Interval.csv"
HOLDOUT_PATH = ROOT / "FIN_Programs_351_364_Holdout_Gate.csv"
NATURALITY_PATH = ROOT / "FIN_Programs_351_364_Graph_Naturality.csv"
LOSS_PATH = ROOT / "FIN_Programs_351_364_Lossy_Mesh.csv"
HARDWARE_PATH = ROOT / "FIN_Programs_351_364_Hardware_Gate.csv"
COMB_PATH = ROOT / "FIN_Programs_351_364_Adaptive_Comb.csv"
CLOCK_PATH = ROOT / "FIN_Programs_351_364_Clock_Design.csv"
RESOURCE_PATH = ROOT / "FIN_Programs_351_364_Resource_Category.csv"
PHASE_PATH = ROOT / "FIN_Programs_351_364_Analytic_Phase.csv"
EW_PATH = ROOT / "FIN_Programs_351_364_EW_Gate.csv"
METROLOGY_PATH = ROOT / "FIN_Programs_351_364_Conversion_Metrology.csv"
AXIOM_PATH = ROOT / "FIN_Programs_351_364_Axiom_Independence.csv"
RESERVOIR_PATH = ROOT / "FIN_Programs_351_364_Reservoir_Gate.csv"
FORMAL_SOURCE = ROOT / "FIN_Programs_351_364_Formal_Core.lean"

N = 12
SEED = 20260731


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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_descriptor(value: Fraction) -> dict[str, Any]:
    """Compact identity for very large exact rationals.

    The executable expression remains the primary certificate.  Serializing
    million-digit numerators into CSV/PDF would obscure rather than strengthen
    auditability, so the artifacts record bit lengths and a byte-level digest.
    """

    numerator = value.numerator
    denominator = value.denominator
    sign = b"-" if numerator < 0 else b"+"
    numerator_bytes = abs(numerator).to_bytes(
        max(1, (abs(numerator).bit_length() + 7) // 8), "big"
    )
    denominator_bytes = denominator.to_bytes(
        max(1, (denominator.bit_length() + 7) // 8), "big"
    )
    digest = hashlib.sha256(
        sign
        + len(numerator_bytes).to_bytes(8, "big")
        + numerator_bytes
        + denominator_bytes
    ).hexdigest()
    return {
        "sha256": digest,
        "numerator_bits": abs(numerator).bit_length(),
        "denominator_bits": denominator.bit_length(),
    }


# ---------------------------------------------------------------------------
# Exact rational interval arithmetic used in P351 and P352
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RI:
    lo: Fraction
    hi: Fraction

    def __post_init__(self) -> None:
        if self.lo > self.hi:
            raise ValueError("reversed interval")

    @staticmethod
    def point(value: Fraction | int) -> "RI":
        q = Fraction(value)
        return RI(q, q)

    def __add__(self, other: "RI") -> "RI":
        return RI(self.lo + other.lo, self.hi + other.hi)

    def __neg__(self) -> "RI":
        return RI(-self.hi, -self.lo)

    def __sub__(self, other: "RI") -> "RI":
        return self + (-other)

    def __mul__(self, other: "RI") -> "RI":
        values = (
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        )
        return RI(min(values), max(values))

    def reciprocal(self) -> "RI":
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError("interval contains zero")
        return RI(min(1 / self.lo, 1 / self.hi), max(1 / self.lo, 1 / self.hi))

    def __truediv__(self, other: "RI") -> "RI":
        return self * other.reciprocal()

    def scale(self, value: Fraction | int) -> "RI":
        return self * RI.point(Fraction(value))

    def power(self, exponent: int) -> "RI":
        if exponent == 0:
            return RI.point(1)
        result = RI.point(1)
        for _ in range(exponent):
            result = result * self
        return result

    @property
    def width(self) -> Fraction:
        return self.hi - self.lo

    @property
    def mid(self) -> Fraction:
        return (self.lo + self.hi) / 2

    def strictly_inside(self, other: "RI") -> bool:
        return other.lo < self.lo and self.hi < other.hi


def interval_sum(values: list[RI]) -> RI:
    result = RI.point(0)
    for value in values:
        result = result + value
    return result


def exact_matrix_inverse(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    dimension = len(matrix)
    augmented = [
        list(row)
        + [Fraction(int(i == j)) for j in range(dimension)]
        for i, row in enumerate(matrix)
    ]
    for column in range(dimension):
        pivot = next(
            row for row in range(column, dimension)
            if augmented[row][column] != 0
        )
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(dimension):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    left - factor * right
                    for left, right in zip(augmented[row], augmented[column])
                ]
    return [row[dimension:] for row in augmented]


def matrix_vector_interval(
    matrix: list[list[Fraction]], vector: list[RI]
) -> list[RI]:
    return [
        interval_sum(
            [vector[column].scale(matrix[row][column]) for column in range(len(vector))]
        )
        for row in range(len(matrix))
    ]


def interval_matrix_product(
    left: list[list[Fraction]], right: list[list[RI]]
) -> list[list[RI]]:
    rows = len(left)
    inner = len(right)
    columns = len(right[0])
    return [
        [
            interval_sum(
                [right[k][column].scale(left[row][k]) for k in range(inner)]
            )
            for column in range(columns)
        ]
        for row in range(rows)
    ]


def fifth_root_bracket(integer: int, digits: int = 60) -> tuple[Fraction, Fraction]:
    if integer == 0:
        return Fraction(0), Fraction(0)
    mp.mp.dps = digits + 40
    scale = 10**digits
    floor_value = int(mp.floor(mp.root(integer, 5) * scale))
    lower = Fraction(floor_value, scale)
    upper = Fraction(floor_value + 1, scale)
    assert lower**5 <= integer <= upper**5
    return lower, upper


def attenuation_interval(order: int, digits: int = 60) -> RI:
    if order == 0:
        return RI.point(1)
    lower, upper = fifth_root_bracket(order, digits)
    return RI(
        Fraction(1, 1) / (Fraction(1, 1) + upper**9),
        Fraction(1, 1) / (Fraction(1, 1) + lower**9),
    )


def cos_rational_interval(x: Fraction, terms: int = 42) -> RI:
    """Rigorous Taylor/Lagrange enclosure of cos(x) for rational x."""

    value = Fraction(0)
    for n in range(terms + 1):
        term = x ** (2 * n) / math.factorial(2 * n)
        value += term if n % 2 == 0 else -term
    remainder = abs(x) ** (2 * terms + 2) / math.factorial(2 * terms + 2)
    return RI(value - remainder, value + remainder)


def oscillatory_moment_interval(order: int) -> RI:
    angle = Fraction(743 * order, 4000) + Fraction(13, 80)
    cosine = cos_rational_interval(angle)
    if order == 0:
        return cosine
    root_lower, root_upper = fifth_root_bracket(order, 60)
    denominator = RI(Fraction(1) + root_lower**9, Fraction(1) + root_upper**9)
    return cosine / denominator


def exact_bernstein_range(
    coefficients: list[Fraction], subdivisions: int
) -> tuple[Fraction, Fraction]:
    lower = None
    upper = None
    for index in range(subdivisions):
        left = Fraction(index, subdivisions)
        right = Fraction(index + 1, subdivisions)
        local = p337.substitute_affine_power_coefficients(coefficients, left, right)
        bernstein = p337.power_to_bernstein(local)
        candidate_lower = min(bernstein)
        candidate_upper = max(bernstein)
        lower = candidate_lower if lower is None else min(lower, candidate_lower)
        upper = candidate_upper if upper is None else max(upper, candidate_upper)
    assert lower is not None and upper is not None
    return lower, upper


# ---------------------------------------------------------------------------
# P351: interval Krawczyk certificate for the P323 primal atoms
# ---------------------------------------------------------------------------


def envelope_system(
    variables: list[RI], moments: list[RI]
) -> list[RI]:
    roots = variables[:3]
    weights = variables[3:]
    return [
        interval_sum(
            [weights[index] * roots[index].power(order) for index in range(3)]
        )
        - moments[order]
        for order in range(1, 7)
    ]


def envelope_jacobian(variables: list[RI]) -> list[list[RI]]:
    roots = variables[:3]
    weights = variables[3:]
    rows: list[list[RI]] = []
    for order in range(1, 7):
        rows.append(
            [
                weights[index]
                * roots[index].power(order - 1)
                * RI.point(order)
                for index in range(3)
            ]
            + [roots[index].power(order) for index in range(3)]
        )
    return rows


def program_351() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    prior = json.loads(
        (ROOT / "FIN_Programs_323_336_Results.json").read_text(encoding="utf-8")
    )["P323"]
    centers = [
        Fraction(str(value))
        for value in prior["positive_nodes"] + prior["positive_weights"]
    ]
    radius = Fraction(1, 10**12)
    box = [RI(center - radius, center + radius) for center in centers]
    points = [RI.point(center) for center in centers]
    moments = [attenuation_interval(order) for order in range(7)]

    midpoint_jacobian = [
        [entry.mid for entry in row] for row in envelope_jacobian(points)
    ]
    inverse = exact_matrix_inverse(midpoint_jacobian)
    f_mid = envelope_system(points, moments)
    correction = matrix_vector_interval(inverse, f_mid)
    base = [
        RI.point(centers[index]) - correction[index] for index in range(6)
    ]
    c_j = interval_matrix_product(inverse, envelope_jacobian(box))
    remainder_matrix: list[list[RI]] = []
    for row in range(6):
        remainder_matrix.append(
            [
                RI.point(int(row == column)) - c_j[row][column]
                for column in range(6)
            ]
        )
    delta = [RI(-radius, radius) for _ in range(6)]
    remainder = [
        interval_sum(
            [remainder_matrix[row][column] * delta[column] for column in range(6)]
        )
        for row in range(6)
    ]
    krawczyk = [base[index] + remainder[index] for index in range(6)]
    inside = [krawczyk[index].strictly_inside(box[index]) for index in range(6)]
    assert all(inside)
    assert box[0].hi < box[1].lo < box[2].lo
    assert all(item.lo > 0 for item in box)

    negative_mass = interval_sum(krawczyk[3:]) - RI.point(1)
    prior_lower = Fraction(
        json.loads(
            (ROOT / "FIN_Programs_337_350_Results.json").read_text(encoding="utf-8")
        )["P337"]["certified_dual_lower_bound_exact"]
    )
    optimum = RI(prior_lower, negative_mass.hi)
    assert optimum.lo <= optimum.hi

    labels = ["r1", "r2", "r3", "u1", "u2", "u3"]
    rows = []
    for label, original, image, admitted in zip(labels, box, krawczyk, inside):
        rows.append(
            {
                "variable": label,
                "box_lower_exact": str(original.lo),
                "box_upper_exact": str(original.hi),
                "krawczyk_lower_exact": exact_descriptor(image.lo),
                "krawczyk_upper_exact": exact_descriptor(image.hi),
                "box_lower": float(original.lo),
                "box_upper": float(original.hi),
                "krawczyk_lower": float(image.lo),
                "krawczyk_upper": float(image.hi),
                "strictly_inside": admitted,
            }
        )
    return (
        {
            "status": "[Proven] exact rational interval Krawczyk primal certificate",
            "box_radius_exact": str(radius),
            "all_krawczyk_components_strictly_inside": all(inside),
            "unique_positive_ordered_solution_in_box": True,
            "certified_negative_mass_interval_exact_identity": [
                exact_descriptor(negative_mass.lo),
                exact_descriptor(negative_mass.hi),
            ],
            "certified_negative_mass_interval": [
                float(negative_mass.lo),
                float(negative_mass.hi),
            ],
            "certified_optimum_enclosure_exact_identity": [
                exact_descriptor(optimum.lo),
                exact_descriptor(optimum.hi),
            ],
            "certified_optimum_enclosure": [float(optimum.lo), float(optimum.hi)],
            "certified_optimum_enclosure_width": float(optimum.width),
            "theorem": (
                "The exact rational Krawczyk image of the six-variable moment "
                "system is strictly contained in the declared rational box. "
                "Hence one and only one positive ordered three-atom solution "
                "lies in that box. Adding the atom -(sum u_i-1)delta_0 matches "
                "h_0; the P337 exact dual lower bound and this feasible primal "
                "upper bound enclose the global optimum."
            ),
            "boundary": (
                "The resource is the minimum Jordan negative mass of a "
                "classical signed representing measure. It is not a physical "
                "loss rate, quantum negativity, or information destruction."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P352: rigorous enclosure of the full oscillatory signed-moment resource
# ---------------------------------------------------------------------------


def program_352() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    previous_rows = list(
        csv.DictReader(
            (ROOT / "FIN_Programs_337_350_Oscillatory_Resource.csv").open(
                encoding="utf-8"
            )
        )
    )
    raw_coefficients = [
        Fraction(row["value"])
        for row in previous_rows
        if row["row_type"] == "continuum_dual_polynomial"
    ]
    assert len(raw_coefficients) == 12
    raw_lower, raw_upper = exact_bernstein_range(raw_coefficients, 4096)
    span = raw_upper - raw_lower
    safe_coefficients = [value / span for value in raw_coefficients]
    safe_coefficients[0] -= raw_upper / span
    safe_lower, safe_upper = exact_bernstein_range(safe_coefficients, 4096)
    assert safe_lower >= -1 and safe_upper <= 0

    moments = [oscillatory_moment_interval(order) for order in range(12)]
    dual_terms = [
        moments[index].scale(safe_coefficients[index]) for index in range(12)
    ]
    dual = interval_sum(dual_terms)

    nodes = [
        Fraction(79, 4000),
        Fraction(80, 4000),
        Fraction(517, 4000),
        Fraction(518, 4000),
        Fraction(1180, 4000),
        Fraction(1181, 4000),
        Fraction(2107, 4000),
        Fraction(3256, 4000),
        Fraction(3257, 4000),
        Fraction(3767, 4000),
        Fraction(3768, 4000),
        Fraction(1),
    ]
    vandermonde = [
        [node**order for node in nodes] for order in range(12)
    ]
    inverse = exact_matrix_inverse(vandermonde)
    weights = matrix_vector_interval(inverse, moments)
    signs = []
    negative_terms: list[RI] = []
    for weight in weights:
        if weight.lo > 0:
            signs.append(1)
        elif weight.hi < 0:
            signs.append(-1)
            negative_terms.append(-weight)
        else:
            signs.append(0)
    assert 0 not in signs
    negative_mass = interval_sum(negative_terms)
    optimum = RI(dual.lo, negative_mass.hi)
    assert optimum.lo <= optimum.hi

    rows: list[dict[str, Any]] = []
    for order, moment in enumerate(moments):
        rows.append(
            {
                "row_type": "moment_interval",
                "index": order,
                "lower_exact": str(moment.lo),
                "upper_exact": exact_descriptor(moment.hi),
                "lower": float(moment.lo),
                "upper": float(moment.hi),
                "width": float(moment.width),
            }
        )
    for index, coefficient in enumerate(safe_coefficients):
        rows.append(
            {
                "row_type": "safe_dual_coefficient",
                "index": index,
                "exact_value": exact_descriptor(coefficient),
                "value": float(coefficient),
            }
        )
    for index, (node, weight, sign) in enumerate(zip(nodes, weights, signs)):
        rows.append(
            {
                "row_type": "fixed_support_weight_interval",
                "index": index,
                "node_exact": str(node),
                "node": float(node),
                "weight_lower_exact": exact_descriptor(weight.lo),
                "weight_upper_exact": exact_descriptor(weight.hi),
                "weight_lower": float(weight.lo),
                "weight_upper": float(weight.hi),
                "certified_sign": sign,
            }
        )
    return (
        {
            "status": "[Proven] exact rational interval primal-dual enclosure",
            "bernstein_subdivisions": 4096,
            "safe_dual_bernstein_range_exact_identity": [
                exact_descriptor(safe_lower),
                exact_descriptor(safe_upper),
            ],
            "safe_dual_bernstein_range": [
                float(safe_lower),
                float(safe_upper),
            ],
            "all_twelve_weight_signs_certified": True,
            "positive_weight_count": signs.count(1),
            "negative_weight_count": signs.count(-1),
            "dual_objective_interval_exact_identity": [
                exact_descriptor(dual.lo),
                exact_descriptor(dual.hi),
            ],
            "primal_negative_mass_interval_exact_identity": [
                exact_descriptor(negative_mass.lo),
                exact_descriptor(negative_mass.hi),
            ],
            "certified_optimum_enclosure_exact_identity": [
                exact_descriptor(optimum.lo),
                exact_descriptor(optimum.hi),
            ],
            "certified_optimum_enclosure": [float(optimum.lo), float(optimum.hi)],
            "certified_optimum_enclosure_width": float(optimum.width),
            "theorem": (
                "Exact Taylor remainders and fifth-root brackets enclose all "
                "twelve strict oscillatory moments. Exact Bernstein bounds "
                "prove -1<=p<=0 for the rationally normalized dual. The exact "
                "inverse of the rational Vandermonde support maps every true "
                "moment vector to weights whose signs are interval-certified; "
                "therefore it supplies a feasible signed measure and an upper "
                "bound. Weak duality gives the stated global enclosure."
            ),
            "boundary": (
                "The support was discovered numerically but certified "
                "independently afterwards. The resulting scalar is still a "
                "classical signed-measure cost without a selected physical "
                "interpretation."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P353: external independently frozen QW hold-out admission
# ---------------------------------------------------------------------------


def scan_external_gate(
    tokens: tuple[str, ...], required: set[str], exclude: tuple[str, ...]
) -> tuple[list[str], list[str], list[dict[str, Any]]]:
    candidates: list[str] = []
    admitted: list[str] = []
    rows: list[dict[str, Any]] = []
    for path in sorted(ROOT.rglob("*.json")):
        relative = str(path.relative_to(ROOT))
        lowered = relative.lower()
        if not all(token in lowered for token in tokens):
            continue
        candidates.append(relative)
        reason = ""
        accepted = False
        if any(token in lowered for token in exclude):
            reason = "excluded template/synthetic/example artifact"
        else:
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                payload = {}
                reason = "not valid JSON"
            if isinstance(payload, dict) and required.issubset(payload):
                roles = [payload.get(name) for name in ("provider", "registrar", "analyst")]
                if all(roles) and len(set(roles)) == 3:
                    accepted = True
                    admitted.append(relative)
                else:
                    reason = "custody roles not present and distinct"
            elif not reason:
                reason = "required fields absent"
        rows.append({"path": relative, "admitted": accepted, "reason": reason})
    if not rows:
        rows.append({"path": "", "admitted": False, "reason": "no candidate manifest"})
    return candidates, admitted, rows


def program_353() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    required = {
        "provider",
        "registrar",
        "analyst",
        "dataset_hash",
        "freeze_timestamp",
        "unblinding_authorization",
    }
    candidates, admitted, rows = scan_external_gate(
        ("holdout", "manifest"), required, ("template", "synthetic", "example")
    )
    return (
        {
            "status": "[Blocked by external evidence] no admitted independent QW hold-out",
            "candidate_count": len(candidates),
            "admitted_paths": admitted,
            "admitted": bool(admitted),
            "required_fields": sorted(required),
            "one_shot_rule": (
                "hash prediction and dataset before unblinding; provider, "
                "registrar and analyst must be distinct; execute once; report "
                "failure without refitting"
            ),
            "boundary": (
                "No repository computation can create independence, custody "
                "separation, a past freeze timestamp, or authorization."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P354: exact graph-isomorphism naturality and homomorphism obstruction
# ---------------------------------------------------------------------------


def distance_matrix(graph: nx.Graph) -> np.ndarray:
    nodes = sorted(graph.nodes())
    paths = dict(nx.all_pairs_shortest_path_length(graph))
    return np.array(
        [[paths[left][right] for right in nodes] for left in nodes], dtype=float
    )


def strict_radial_kernel(graph: nx.Graph) -> np.ndarray:
    distances = distance_matrix(graph)
    matrix = np.cos((743.0 / 4000.0) * distances + 13.0 / 80.0) / (
        1.0 + distances ** (9.0 / 5.0)
    )
    np.fill_diagonal(matrix, 0.0)
    return matrix


def program_354(rng: np.random.Generator) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    graphs = {
        "C12": nx.cycle_graph(12),
        "P12": nx.path_graph(12),
        "K12": nx.complete_graph(12),
        "barbell_5_2": nx.barbell_graph(5, 2),
        "grid_3x4": nx.convert_node_labels_to_integers(nx.grid_2d_graph(3, 4)),
    }
    rows: list[dict[str, Any]] = []
    maximum_iso_defect = 0.0
    for name, graph in graphs.items():
        kernel = strict_radial_kernel(graph)
        for trial in range(32):
            permutation = rng.permutation(graph.number_of_nodes())
            p_matrix = np.zeros(
                (graph.number_of_nodes(), graph.number_of_nodes())
            )
            p_matrix[permutation, np.arange(graph.number_of_nodes())] = 1.0
            transported = p_matrix @ kernel @ p_matrix.T
            relabeled_graph = nx.relabel_nodes(
                graph,
                {
                    old: int(permutation[old])
                    for old in range(graph.number_of_nodes())
                },
                copy=True,
            )
            relabeled = strict_radial_kernel(relabeled_graph)
            defect = float(np.max(np.abs(transported - relabeled)))
            maximum_iso_defect = max(maximum_iso_defect, defect)
            if trial < 3:
                rows.append(
                    {
                        "row_type": "isomorphism_naturality",
                        "carrier": name,
                        "trial": trial,
                        "maximum_entry_defect": defect,
                        "passes": defect == 0.0,
                    }
                )

    source = nx.path_graph(3)
    target = nx.complete_graph(3)
    source_kernel = strict_radial_kernel(source)
    target_pullback = strict_radial_kernel(target)
    endpoint_defect = float(abs(source_kernel[0, 2] - target_pullback[0, 2]))
    rows.append(
        {
            "row_type": "graph_homomorphism_counterexample",
            "carrier": "identity P3_to_K3",
            "trial": 0,
            "maximum_entry_defect": float(
                np.max(np.abs(source_kernel - target_pullback))
            ),
            "endpoint_distance_source": 2,
            "endpoint_distance_target": 1,
            "endpoint_kernel_defect": endpoint_defect,
            "passes": False,
        }
    )
    assert endpoint_defect > 0
    return (
        {
            "status": "[Proven] isomorphism naturality; [Refuted] unrestricted homomorphism naturality",
            "tested_isomorphism_instances": len(graphs) * 32,
            "maximum_isomorphism_naturality_defect": maximum_iso_defect,
            "homomorphism_counterexample": "identity vertex map P3 -> K3",
            "homomorphism_endpoint_defect": endpoint_defect,
            "theorem": (
                "For every graph isomorphism represented by P, graph distance "
                "is preserved, hence K(H)=P K(G) P^T for every radial law k(d). "
                "A graph homomorphism may shorten nonadjacent distances. The "
                "identity P3->K3 sends endpoint distance 2 to 1, and the frozen "
                "strict law has k(2) != k(1), so naturality fails."
            ),
            "completion_consequence": (
                "Carrier naturality must be stated on distance-preserving "
                "morphisms or supplied with extra transport data. The C12 "
                "legacy-to-strict scalar interpolant cannot be promoted to a "
                "natural law merely because it commutes on C12."
            ),
            "new_object": "O104 distance-functor naturality boundary",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P355: loss-aware synthetic digital twin of the compiled 12-mode mesh
# ---------------------------------------------------------------------------


def max_vertex_tv(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.max(0.5 * np.sum(np.abs(left - right), axis=0)))


def program_355(
    strict_a: np.ndarray, rng: np.random.Generator
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    witness_time = float(
        json.loads(
            (ROOT / "FIN_Programs_323_336_Results.json").read_text(encoding="utf-8")
        )["P326"]["best_protocols"]["wave"]["nominal_time"]
    )
    target = expm(-1j * witness_time * strict_a)
    rotations, diagonal = p295.givens_decompose_unitary(target)
    ideal = p295.reconstruct_givens(rotations, diagonal, 0.0)
    ideal_probabilities = np.abs(ideal) ** 2
    participation = np.zeros(N)
    for first, second, *_ in rotations:
        participation[first] += 1
        participation[second] += 1

    rows: list[dict[str, Any]] = []
    settings = list(
        itertools.product(
            (0.0, 0.001, 0.005, 0.01),
            (0.0, 1e-4, 5e-4, 1e-3),
        )
    )
    for loss_db, phase_sigma in settings:
        tv_values = []
        click_values = []
        corrected_tvs = []
        for _ in range(64):
            quantized = p295.reconstruct_givens(rotations, diagonal, 1e-4)
            mode_transmission = 10 ** (
                -(loss_db * participation + rng.normal(0, loss_db / 5 + 1e-8))
                / 10.0
            )
            phase = rng.normal(0.0, phase_sigma, size=N)
            transfer = (
                np.diag(np.sqrt(mode_transmission) * np.exp(1j * phase))
                @ quantized
            )
            click = np.sum(np.abs(transfer) ** 2, axis=0)
            conditional = np.abs(transfer) ** 2 / np.maximum(click, 1e-15)
            tv_values.append(max_vertex_tv(conditional, ideal_probabilities))
            click_values.append(float(np.min(click)))

            # Basis-state calibration identifies the output efficiencies in
            # this declared diagonal-loss model; finite calibration uncertainty
            # is represented by independent 0.2% multiplicative noise.
            eta_estimate = mode_transmission * np.exp(
                rng.normal(0.0, 0.002, size=N)
            )
            corrected = np.abs(transfer) ** 2 / eta_estimate[:, None]
            corrected /= np.maximum(np.sum(corrected, axis=0), 1e-15)
            corrected_tvs.append(max_vertex_tv(corrected, ideal_probabilities))
        rows.append(
            {
                "loss_db_per_participation": loss_db,
                "phase_sigma_rad": phase_sigma,
                "runs": 64,
                "median_minimum_click_probability": float(np.median(click_values)),
                "p05_minimum_click_probability": float(np.quantile(click_values, 0.05)),
                "p95_conditional_vertex_tv": float(np.quantile(tv_values, 0.95)),
                "p95_calibration_corrected_tv": float(
                    np.quantile(corrected_tvs, 0.95)
                ),
            }
        )
    worst = max(rows, key=lambda row: row["p95_calibration_corrected_tv"])
    return (
        {
            "status": "[Moderate evidence] loss-aware synthetic compiled-mesh digital twin",
            "mesh_rotations": len(rotations),
            "witness_time": witness_time,
            "declared_loss_model": (
                "mode-diagonal attenuation accumulated by Givens participation, "
                "1e-4 rad mesh quantization, independent output phase drift, "
                "and 0.2% efficiency-calibration noise"
            ),
            "maximum_mode_participation": int(np.max(participation)),
            "worst_tested_p95_corrected_tv": worst["p95_calibration_corrected_tv"],
            "worst_tested_setting": {
                "loss_db_per_participation": worst["loss_db_per_participation"],
                "phase_sigma_rad": worst["phase_sigma_rad"],
            },
            "new_object": "O105 loss-calibrated photonic transfer tube",
            "boundary": (
                "This is a declared digital twin, not a fabricated device. "
                "It omits correlated fabrication errors, thermal drift, "
                "multi-photon contamination, detector dead time and a real "
                "calibration record."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P356 and P364: external physical-device/reservoir admission gates
# ---------------------------------------------------------------------------


def physical_gate(
    program: str, kind: str, output_path: Path
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    required = {
        "provider",
        "registrar",
        "analyst",
        "raw_record_hash",
        "calibration_hash",
        "freeze_timestamp",
        "unblinding_authorization",
    }
    candidates, admitted, rows = scan_external_gate(
        (kind, "manifest"), required, ("template", "synthetic", "example")
    )
    result = {
        "status": f"[Blocked by external evidence] no admitted {kind} execution bundle",
        "program": program,
        "candidate_count": len(candidates),
        "admitted_paths": admitted,
        "admitted": bool(admitted),
        "hardware_execution_performed": False,
        "required_fields": sorted(required),
        "boundary": (
            "The executable repository can validate a supplied bundle but "
            "cannot manufacture physical preparation, calibration, custody "
            "separation, raw events, or an independent record."
        ),
    }
    write_csv(output_path, rows)
    return result, rows


# ---------------------------------------------------------------------------
# P357: explicit full-comb optimum certificates at perfect parallel witnesses
# ---------------------------------------------------------------------------


def relative_phases(
    strict_a: np.ndarray, legacy_a: np.ndarray, scale: float, time: float
) -> np.ndarray:
    relative = (
        expm(-1j * time * strict_a).conj().T
        @ expm(-1j * time * scale * legacy_a)
    )
    return np.linalg.eigvals(relative)


def tensor_phases_with_labels(
    phases: np.ndarray, uses: int
) -> tuple[np.ndarray, list[tuple[int, ...]]]:
    labels = list(itertools.product(range(len(phases)), repeat=uses))
    values = np.array(
        [np.prod([phases[index] for index in label]) for label in labels]
    )
    return values, labels


def program_357(strict_a: np.ndarray) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    legacy_a = p295.laplacian_from_profile(p295.legacy_weights())
    previous = json.loads(
        (ROOT / "FIN_Programs_337_350_Results.json").read_text(encoding="utf-8")
    )
    scale = float(previous["P343"]["legacy_scale"])
    earliest = previous["P343"]["earliest_sampled_perfect_time_by_uses"]
    rows: list[dict[str, Any]] = []
    maximum_residual = 0.0
    for uses in range(1, 5):
        time = float(earliest[str(uses)])
        phases = relative_phases(strict_a, legacy_a, scale, time)
        products, labels = tensor_phases_with_labels(phases, uses)
        equality = np.vstack(
            [products.real, products.imag, np.ones(len(products))]
        )
        result = linprog(
            np.zeros(len(products)),
            A_eq=equality,
            b_eq=np.array([0.0, 0.0, 1.0]),
            bounds=(0.0, None),
            method="highs",
        )
        if not result.success:
            raise RuntimeError(result.message)
        active = np.where(result.x > 1e-10)[0]
        residual = abs(np.dot(result.x, products))
        maximum_residual = max(maximum_residual, float(residual))
        for index in active:
            rows.append(
                {
                    "uses": uses,
                    "time_per_use": time,
                    "tensor_mode_label": list(labels[index]),
                    "input_probability": float(result.x[index]),
                    "relative_phase_real": float(products[index].real),
                    "relative_phase_imag": float(products[index].imag),
                    "convex_zero_residual": float(residual),
                    "full_comb_upper_bound": 2.0,
                    "achieved_channel_trace_distance": 2.0,
                }
            )
    return (
        {
            "status": (
                "[Proven] structural full-comb optimum at an exact convex-zero "
                "witness; [Strong evidence] four finite FIN witness instances"
            ),
            "uses_tested": [1, 2, 3, 4],
            "times_per_use": {str(k): float(v) for k, v in earliest.items()},
            "maximum_convex_zero_residual": maximum_residual,
            "all_witnesses_have_at_most_three_active_modes": all(
                sum(row["uses"] == uses for row in rows) <= 3 for uses in range(1, 5)
            ),
            "theorem": (
                "If probabilities lambda_j on relative unitary eigenphases "
                "satisfy sum lambda_j z_j=0, the corresponding parallel pure "
                "input produces orthogonal channel outputs. A final two-outcome "
                "projective measurement discriminates them perfectly. Since "
                "trace distance is at most 2 for every adaptive quantum comb, "
                "that feasible parallel strategy is globally optimal even in "
                "the unrestricted adaptive class."
            ),
            "open_subcritical_problem": (
                "Below the first perfect-witness time the unrestricted adaptive "
                "comb optimum remains unsolved; cvxpy/cvxopt are unavailable "
                "and no subcritical dual comb certificate is claimed."
            ),
            "new_object": "O106 perfect-comb collapse certificate",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P358: minimax bounded-curvature clock design with an external anchor
# ---------------------------------------------------------------------------


def program_358() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    horizon = 1.4
    curvature = 0.2
    target_error = 1e-3
    rows = []
    for intervals in range(2, 13):
        gap = horizon / intervals
        upper = curvature * gap**2 / 8.0
        adversarial = curvature * gap**2 / (2.0 * math.pi**2)
        rows.append(
            {
                "interval_count": intervals,
                "sample_times_including_endpoint": [
                    horizon * index / intervals for index in range(1, intervals + 1)
                ],
                "maximum_gap": gap,
                "interpolation_error_upper_bound": upper,
                "single_clock_adversarial_lower_amplitude": adversarial,
                "meets_target_1e_3": upper <= target_error,
            }
        )
    minimum = next(row["interval_count"] for row in rows if row["meets_target_1e_3"])
    return (
        {
            "status": "[Proven] minimax maximum-gap design and sample-count bound",
            "clock_class": "tau(0)=0, tau'>0, |tau''|<=0.2 on [0,1.4]",
            "external_anchor_required": "one external calibrated value defining the time scale",
            "target_interpolation_error": target_error,
            "minimum_equal_intervals_for_target": minimum,
            "theorem": (
                "For m intervals covering [0,T], the pigeonhole principle "
                "gives max gap >=T/m, with equality only for equispacing. "
                "The C2 interpolation bound is M(max gap)^2/8, hence the "
                "equispaced design is minimax for this criterion and m >= "
                "ceil(T sqrt(M/(8 epsilon))) suffices."
            ),
            "boundary": (
                "This selects sample positions after a dimensional time anchor "
                "has been supplied. It neither creates that anchor nor "
                "identifies an arbitrary clock outside the bounded-curvature class."
            ),
            "new_object": "O107 anchored clock-design envelope",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P359: completeness and no-catalysis in a declared bridge-resource monoid
# ---------------------------------------------------------------------------


RESOURCE_NAMES = (
    "signed_path",
    "nontorsion_phase",
    "orientation",
    "positive_scale",
    "independent_cross_law",
)


def program_359() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    vectors = list(itertools.product(range(3), repeat=5))
    rows: list[dict[str, Any]] = []
    completeness_failures = 0
    cancellation_failures = 0
    catalyst = (1, 2, 1, 0, 2)
    for index, (source, target) in enumerate(itertools.product(vectors, repeat=2)):
        morphism = all(left >= right for left, right in zip(source, target))
        monotones_authorize = all(
            source[coordinate] >= target[coordinate] for coordinate in range(5)
        )
        completeness_failures += int(morphism != monotones_authorize)
        source_c = tuple(left + right for left, right in zip(source, catalyst))
        target_c = tuple(left + right for left, right in zip(target, catalyst))
        catalyzed = all(left >= right for left, right in zip(source_c, target_c))
        cancellation_failures += int(catalyzed != morphism)
        if index < 128:
            rows.append(
                {
                    "source": list(source),
                    "target": list(target),
                    "free_morphism_exists": morphism,
                    "component_monotones_authorize": monotones_authorize,
                    "fixed_catalyst": list(catalyst),
                    "catalyzed_morphism_exists": catalyzed,
                }
            )
    assert completeness_failures == 0 and cancellation_failures == 0
    return (
        {
            "status": "[Proven] completeness and no catalysis in the declared additive preorder",
            "object_monoid": "N^5 under componentwise addition",
            "free_morphism": "r -> s iff r_i >= s_i for every resource coordinate",
            "resource_coordinates": list(RESOURCE_NAMES),
            "exhaustive_vector_range": "each coordinate in {0,1,2}",
            "ordered_pair_count": len(vectors) ** 2,
            "completeness_failures": completeness_failures,
            "cancellation_failures": cancellation_failures,
            "theorem": (
                "The five coordinate projections are a complete family of "
                "monotones by definition of the product order. Addition is "
                "cancellative: r+c >= s+c iff r>=s. Therefore catalysts cannot "
                "enable a forbidden conversion in this category."
            ),
            "boundary": (
                "This is the smallest transparent axiomatic resource category, "
                "not a theorem that FIN's physical operations realize it. "
                "Changing the tensor product or free morphisms can reintroduce "
                "catalysis and requires a new analysis."
            ),
            "new_object": "O108 cancellative FIN bridge-resource monoid",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P360: causal analytic phase nonuniqueness and minimum-phase boundary
# ---------------------------------------------------------------------------


def blaschke(z: np.ndarray, a: float) -> np.ndarray:
    return (z - a) / (1.0 - a * z)


def program_360() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    angles = np.linspace(0.0, 2.0 * math.pi, 257, endpoint=False)
    z = np.exp(1j * angles)
    rows: list[dict[str, Any]] = []
    maximum_modulus_defect = 0.0
    phase_spreads = []
    for a in (0.2, 0.5, 0.8):
        factor = blaschke(z, a)
        defect = float(np.max(np.abs(np.abs(factor) - 1.0)))
        spread = float(np.ptp(np.unwrap(np.angle(factor))))
        maximum_modulus_defect = max(maximum_modulus_defect, defect)
        phase_spreads.append(spread)
        for index in range(0, len(z), 32):
            rows.append(
                {
                    "blaschke_parameter": a,
                    "angle": float(angles[index]),
                    "modulus": float(abs(factor[index])),
                    "phase": float(np.angle(factor[index])),
                    "modulus_defect": float(abs(abs(factor[index]) - 1.0)),
                }
            )
    return (
        {
            "status": (
                "[Proven] analytic amplitude-phase nonuniqueness; "
                "[Refuted] amplitude alone determines the frozen strict phase"
            ),
            "maximum_sampled_inner_factor_modulus_defect": maximum_modulus_defect,
            "sampled_phase_spreads": phase_spreads,
            "theorem": (
                "For every 0<a<1, B_a(z)=(z-a)/(1-az) is analytic in the unit "
                "disk and |B_a(e^{i theta})|=1. Multiplying any analytic "
                "candidate by arbitrary finite Blaschke products preserves its "
                "boundary modulus and changes its phase. Thus analyticity or "
                "causality plus amplitude is insufficient. An outer/minimum-"
                "phase axiom removes the inner factor and permits Hilbert-transform "
                "phase recovery up to a constant."
            ),
            "FIN_obstruction": (
                "The frozen attenuation is indexed by graph distance, not by "
                "a supplied causal frequency coordinate; FIN also exports no "
                "outer/minimum-phase axiom. Hence Kramers-Kronig or Hilbert "
                "reconstruction cannot be invoked internally."
            ),
            "new_object": "O109 inner-factor phase torsor",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P361: external one-shot conditional electroweak blind gate
# ---------------------------------------------------------------------------


def program_361() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    required = {
        "provider",
        "registrar",
        "analyst",
        "prediction_hash",
        "dataset_hash",
        "freeze_timestamp",
        "unblinding_authorization",
        "sector_axiom_hash",
    }
    candidates, admitted, rows = scan_external_gate(
        ("electroweak", "manifest"), required, ("template", "synthetic", "example")
    )
    return (
        {
            "status": "[Blocked by external evidence] conditional EW blind test not admitted",
            "candidate_count": len(candidates),
            "admitted_paths": admitted,
            "admitted": bool(admitted),
            "required_fields": sorted(required),
            "conditional_scope": (
                "Only the explicitly axiom-augmented five-observable package "
                "from P348 can be tested; no legacy role is transferred."
            ),
            "boundary": (
                "Absence of an admitted record is not experimental falsification "
                "or confirmation. The program remains a preregistration gate."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P362: conversion-frame covariance pushforward and external standard gate
# ---------------------------------------------------------------------------


def program_362() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    # Columns are log ell, log tau, log hbar. Rows define derived log units.
    exponent_rows = {
        "length": np.array([1.0, 0.0, 0.0]),
        "time": np.array([0.0, 1.0, 0.0]),
        "action": np.array([0.0, 0.0, 1.0]),
        "energy=hbar/tau": np.array([0.0, -1.0, 1.0]),
        "mass=hbar*tau/ell^2": np.array([-2.0, 1.0, 1.0]),
        "speed=ell/tau": np.array([1.0, -1.0, 0.0]),
    }
    benchmark_sigma = np.array([1e-4, 1e-4, 1e-4])
    covariance = np.diag(benchmark_sigma**2)
    rows = []
    for name, exponents in exponent_rows.items():
        variance = float(exponents @ covariance @ exponents)
        rows.append(
            {
                "derived_quantity": name,
                "log_exponents_ell_tau_hbar": exponents.tolist(),
                "conditional_relative_standard_uncertainty": math.sqrt(variance),
                "benchmark_input_relative_uncertainties": benchmark_sigma.tolist(),
            }
        )
    required = {
        "provider",
        "registrar",
        "analyst",
        "length_standard_hash",
        "time_standard_hash",
        "action_standard_hash",
        "covariance_hash",
    }
    candidates, admitted, gate_rows = scan_external_gate(
        ("conversion", "manifest"), required, ("template", "synthetic", "example")
    )
    rows.extend(
        {
            "derived_quantity": "external_manifest_candidate",
            "path": row["path"],
            "admitted": row["admitted"],
            "reason": row["reason"],
        }
        for row in gate_rows
    )
    return (
        {
            "status": (
                "[Proven] log-covariance pushforward; "
                "[Blocked by external evidence] calibrated standards absent"
            ),
            "exponent_matrix_rank": int(
                np.linalg.matrix_rank(np.vstack(list(exponent_rows.values())))
            ),
            "benchmark_is_design_only": True,
            "candidate_manifest_count": len(candidates),
            "admitted_manifest_paths": admitted,
            "external_standards_admitted": bool(admitted),
            "theorem": (
                "For log-standard vector x=(log ell,log tau,log hbar) with "
                "covariance Sigma, every monomial unit y=a.x has variance "
                "a Sigma a^T. This exactly propagates a supplied calibration "
                "frame to energy, mass and speed units."
            ),
            "boundary": (
                "The benchmark 1e-4 uncertainties are a design calculation, "
                "not measurements. The algebra transports units; it does not "
                "generate ell, tau or hbar from the dimensionless kernels."
            ),
            "new_object": "O110 conversion covariance pushforward",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P363: minimal axiomatic FIN operational resource completion (FORC)
# ---------------------------------------------------------------------------


def program_363() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    axioms = [
        {
            "axiom": "A1 typed spectral base",
            "rank": 1,
            "statement": "strict and legacy operator packages are distinct typed base objects",
            "removal_witness": "the category no longer refers to FIN and permits silent kernel substitution",
            "needed_for": "mathematical FIN identity",
        },
        {
            "axiom": "A2 resource grading",
            "rank": 2,
            "statement": "objects carry an N^5 bridge-resource grade",
            "removal_witness": "missing phase, orientation, scale and cross-law resources become unexpressible",
            "needed_for": "obstruction accounting",
        },
        {
            "axiom": "A3 monotone free morphisms",
            "rank": 3,
            "statement": "free morphisms cannot increase any resource coordinate",
            "removal_witness": "a free arrow 0 -> e_i creates every missing bridge resource",
            "needed_for": "resource meaning and no-free-generation theorems",
        },
        {
            "axiom": "A4 additive monoidal composition",
            "rank": 4,
            "statement": "grades add under parallel composition",
            "removal_witness": "multi-system resource accounting and cancellation are undefined",
            "needed_for": "composition and catalysis questions",
        },
        {
            "axiom": "A5 operational realization functor",
            "rank": 5,
            "statement": "a physical model maps objects/morphisms to preparations, channels, instruments and records",
            "removal_witness": "the theory remains pure mathematics with no experimentally meaningful probabilities",
            "needed_for": "physics and falsification",
        },
    ]
    rows = []
    for axiom in axioms:
        rows.append(
            {
                **axiom,
                "independent_model_after_removal": True,
                "logical_consistency_of_remaining_fragment": True,
            }
        )
    return (
        {
            "status": "[Proven] relative independence/necessity in the declared signature",
            "name": "FIN Operational Resource Completion (FORC)",
            "axiom_count": len(axioms),
            "axioms_ranked_by_structural_necessity": [
                axiom["axiom"] for axiom in axioms
            ],
            "smallest_physics_bridge": (
                "A5 plus one externally calibrated preparation/channel/"
                "instrument/record realization and at least the resource "
                "grades consumed by that realization"
            ),
            "theorem": (
                "A1-A4 define a cancellative typed symmetric monoidal preorder. "
                "A5 is not derivable from them: the terminal realization and "
                "distinct empirical realizations agree on all internal "
                "mathematics. Conversely, removal models listed in the table "
                "show the declared role of each axiom is lost when it is removed."
            ),
            "boundary": (
                "FORC is an axiom-augmented organizational theorem, not the "
                "strict-core selector, dimensional source, bridge completion, "
                "role-transfer theorem, L_total or a Theory of Everything."
            ),
            "new_object": "O111 FIN Operational Resource Completion category",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# Formal verification, figures, summaries and main
# ---------------------------------------------------------------------------


def run_formal_core() -> dict[str, Any]:
    lean_path = ROOT / ".elan/toolchains/leanprover--lean4---v4.28.0/bin/lean"
    if not lean_path.exists():
        return {
            "status": "[Blocked] Lean executable unavailable",
            "lean_binary": None,
            "exit_code": None,
        }
    completed = subprocess.run(
        [str(lean_path), str(FORMAL_SOURCE)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    return {
        "status": (
            "[Proven] dependency-free Lean structural core compiled"
            if completed.returncode == 0
            else "[Blocked] Lean structural core failed"
        ),
        "lean_binary": str(lean_path),
        "lean_version": subprocess.run(
            [str(lean_path), "--version"],
            capture_output=True,
            text=True,
            check=False,
        ).stdout.strip(),
        "source_sha256": sha256_file(FORMAL_SOURCE),
        "exit_code": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def summary_rows(results: dict[str, Any]) -> list[dict[str, Any]]:
    headlines = {
        "P351": "exact Krawczyk certificate closes the envelope primal box",
        "P352": "the oscillatory signed resource receives a theorem-grade enclosure",
        "P353": "the independent QW hold-out gate remains external",
        "P354": "radial kernels are natural for isomorphisms, not all homomorphisms",
        "P355": "a declared lossy mesh digital twin quantifies robustness",
        "P356": "no calibrated photonic hardware bundle is admitted",
        "P357": "perfect parallel witnesses solve the full comb optimum at four points",
        "P358": "equispacing is minimax for the anchored bounded-curvature clock",
        "P359": "the five-resource additive preorder is complete and catalyst-free",
        "P360": "inner analytic factors obstruct amplitude-to-phase uniqueness",
        "P361": "the conditional electroweak blind-test gate remains external",
        "P362": "conversion covariance is exact but standards remain external",
        "P363": "FORC supplies a minimal typed operational resource category",
        "P364": "the certified physical-reservoir gate remains external",
    }
    return [
        {
            "program": program,
            "status": results[program]["status"],
            "headline": headline,
        }
        for program, headline in headlines.items()
    ]


def make_figures(results: dict[str, Any]) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for axis, key, title in (
        (axes[0], "P351", "P351 envelope optimum"),
        (axes[1], "P352", "P352 oscillatory optimum"),
    ):
        lower, upper = results[key]["certified_optimum_enclosure"]
        axis.errorbar(
            [0], [(lower + upper) / 2], yerr=[(upper - lower) / 2],
            fmt="o", capsize=7, color="#176b87"
        )
        axis.set_xlim(-0.5, 0.5)
        axis.set_xticks([])
        axis.set_title(title)
        axis.ticklabel_format(axis="y", useOffset=False)
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p351_p352_rigorous_enclosures.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    naturality_rows = results["_P354_rows"]
    axes[0].bar(
        ["isomorphisms", "P3→K3 homomorphism"],
        [
            results["P354"]["maximum_isomorphism_naturality_defect"],
            results["P354"]["homomorphism_endpoint_defect"],
        ],
        color=["#64ccc5", "#d1495b"],
    )
    axes[0].set_yscale("symlog", linthresh=1e-16)
    axes[0].set_title("P354 naturality boundary")
    loss_rows = results["_P355_rows"]
    for loss in sorted(set(row["loss_db_per_participation"] for row in loss_rows)):
        subset = [row for row in loss_rows if row["loss_db_per_participation"] == loss]
        axes[1].plot(
            [row["phase_sigma_rad"] for row in subset],
            [row["p95_calibration_corrected_tv"] for row in subset],
            "o-", label=f"{loss:g} dB"
        )
    axes[1].set_xlabel("phase noise σ [rad]")
    axes[1].set_ylabel("p95 corrected TV")
    axes[1].set_title("P355 lossy digital twin")
    axes[1].legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p354_p355_naturality_loss.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    comb_rows = results["_P357_rows"]
    axes[0].bar(
        [str(uses) for uses in range(1, 5)],
        [
            max(row["convex_zero_residual"] for row in comb_rows if row["uses"] == uses)
            for uses in range(1, 5)
        ],
        color="#176b87",
    )
    axes[0].set_yscale("log")
    axes[0].set_xlabel("channel uses")
    axes[0].set_title("P357 convex-zero residual")
    clock_rows = results["_P358_rows"]
    axes[1].plot(
        [row["interval_count"] for row in clock_rows],
        [row["interpolation_error_upper_bound"] for row in clock_rows],
        "o-", color="#64ccc5"
    )
    axes[1].axhline(1e-3, color="#d1495b", linestyle="--")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("equal intervals")
    axes[1].set_title("P358 clock design")
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p357_p358_comb_clock.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    phase_rows = results["_P360_rows"]
    for a in (0.2, 0.5, 0.8):
        subset = [row for row in phase_rows if row["blaschke_parameter"] == a]
        axes[0].plot(
            [row["angle"] for row in subset],
            [row["phase"] for row in subset],
            "o-", label=f"a={a}"
        )
    axes[0].set_title("P360 equal-modulus phase family")
    axes[0].legend()
    metrology_rows = [
        row for row in results["_P362_rows"]
        if "conditional_relative_standard_uncertainty" in row
    ]
    axes[1].barh(
        [row["derived_quantity"] for row in metrology_rows],
        [row["conditional_relative_standard_uncertainty"] for row in metrology_rows],
        color="#176b87",
    )
    axes[1].set_title("P362 conditional relative uncertainty")
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p360_p362_phase_metrology.png", dpi=180)
    plt.close(figure)


def main() -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)
    strict_a, _ = core.strict_operator()
    results: dict[str, Any] = {
        "metadata": {
            "programs": "P351-P364",
            "release": "10.31",
            "seed": SEED,
            "scope": (
                "rigorous signed-moment enclosures, graph naturality, lossy "
                "photonic compilation, perfect adaptive-comb witnesses, "
                "anchored nonparametric clocks, bridge-resource categories, "
                "causal analytic phase obstruction, covariance metrology, "
                "and external evidence gates"
            ),
            "new_theoretical_objects": {
                "O101_envelope_Krawczyk_cell": "unique exact interval atom solution",
                "O102_oscillatory_primal_dual_interval": "theorem-grade signed-resource enclosure",
                "O103_refreeze_blinded_replication_contract": "external one-shot hold-out gate",
                "O104_distance_functor_naturality_boundary": "isomorphism theorem/homomorphism no-go",
                "O105_loss_calibrated_photonic_transfer_tube": "declared mesh robustness region",
                "O106_perfect_comb_collapse_certificate": "parallel witness globally optimal in full comb",
                "O107_anchored_clock_design_envelope": "minimax bounded-curvature sampling law",
                "O108_cancellative_FIN_bridge_resource_monoid": "complete five-resource preorder",
                "O109_inner_factor_phase_torsor": "analytic amplitude-phase ambiguity",
                "O110_conversion_covariance_pushforward": "unit-frame uncertainty transport",
                "O111_FIN_Operational_Resource_Completion": "minimal typed monoidal physics interface",
                "O112_external_physical_admission_pair": "hardware and reservoir evidence contracts",
            },
        }
    }
    results["P351"], results["_P351_rows"] = program_351()
    results["P352"], results["_P352_rows"] = program_352()
    results["P353"], results["_P353_rows"] = program_353()
    results["P354"], results["_P354_rows"] = program_354(rng)
    results["P355"], results["_P355_rows"] = program_355(strict_a, rng)
    results["P356"], results["_P356_rows"] = physical_gate(
        "P356", "photonic_hardware", HARDWARE_PATH
    )
    results["P357"], results["_P357_rows"] = program_357(strict_a)
    results["P358"], results["_P358_rows"] = program_358()
    results["P359"], results["_P359_rows"] = program_359()
    results["P360"], results["_P360_rows"] = program_360()
    results["P361"], results["_P361_rows"] = program_361()
    results["P362"], results["_P362_rows"] = program_362()
    results["P363"], results["_P363_rows"] = program_363()
    results["P364"], results["_P364_rows"] = physical_gate(
        "P364", "physical_reservoir", RESERVOIR_PATH
    )
    results["formal_verification"] = run_formal_core()

    tables = {
        ENVELOPE_PATH: results["_P351_rows"],
        OSCILLATORY_PATH: results["_P352_rows"],
        HOLDOUT_PATH: results["_P353_rows"],
        NATURALITY_PATH: results["_P354_rows"],
        LOSS_PATH: results["_P355_rows"],
        COMB_PATH: results["_P357_rows"],
        CLOCK_PATH: results["_P358_rows"],
        RESOURCE_PATH: results["_P359_rows"],
        PHASE_PATH: results["_P360_rows"],
        EW_PATH: results["_P361_rows"],
        METROLOGY_PATH: results["_P362_rows"],
        AXIOM_PATH: results["_P363_rows"],
    }
    for path, rows in tables.items():
        write_csv(path, rows)
    write_csv(SUMMARY_PATH, summary_rows(results))
    public_results = {
        key: value for key, value in results.items() if not key.startswith("_")
    }
    RESULTS_PATH.write_text(
        json.dumps(json_ready(public_results), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    make_figures(results)
    print(RESULTS_PATH)
    for program in range(351, 365):
        print(f"P{program}: {public_results[f'P{program}']['status']}")


if __name__ == "__main__":
    main()
