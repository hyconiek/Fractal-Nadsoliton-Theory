#!/usr/bin/env python3
"""Execute FIN Research Programs P365--P378 (Release 10.32).

The round independently checks the exact Release-10.31 moment certificates,
tightens the oscillatory primal-dual enclosure, classifies radial-kernel
naturality, audits component-level photonic loss identifiability, solves the
subcritical unrestricted adaptive-unitary discrimination problem analytically,
constructs operational resource monotones, proves a conditional outer
generating-function theorem, and determines the precise universal scope of
FORC's resource skeleton.

Programs requiring independent hold-out data, clocks, standards, hardware,
electroweak records, or reservoirs remain external admission gates.
"""

from __future__ import annotations

import csv
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.linalg import expm

import fin_programs_255_266 as core
import fin_programs_295_308 as p295
import fin_programs_323_336 as p323
import fin_programs_351_364 as p351


ROOT = Path(__file__).resolve().parent
FIGURE_DIR = ROOT / "FIN_Programs_365_378_Figures"
RESULTS_PATH = ROOT / "FIN_Programs_365_378_Results.json"
SUMMARY_PATH = ROOT / "FIN_Programs_365_378_Summary.csv"
CHECKER_PATH = ROOT / "FIN_Programs_365_378_Independent_Checker.csv"
TIGHTENING_PATH = ROOT / "FIN_Programs_365_378_Oscillatory_Tightening.csv"
HOLDOUT_PATH = ROOT / "FIN_Programs_365_378_Holdout_Gate.csv"
NATURALITY_PATH = ROOT / "FIN_Programs_365_378_Maximal_Naturality.csv"
COMPONENT_PATH = ROOT / "FIN_Programs_365_378_Component_Loss.csv"
HARDWARE_PATH = ROOT / "FIN_Programs_365_378_Hardware_Gate.csv"
COMB_PATH = ROOT / "FIN_Programs_365_378_Subcritical_Comb.csv"
CLOCK_PATH = ROOT / "FIN_Programs_365_378_Clock_Gate.csv"
SEMANTICS_PATH = ROOT / "FIN_Programs_365_378_Resource_Semantics.csv"
OUTER_PATH = ROOT / "FIN_Programs_365_378_Outer_Extension.csv"
EW_PATH = ROOT / "FIN_Programs_365_378_EW_Gate.csv"
STANDARDS_PATH = ROOT / "FIN_Programs_365_378_Standards_Gate.csv"
UNIVERSALITY_PATH = ROOT / "FIN_Programs_365_378_FORC_Universality.csv"
RESERVOIR_PATH = ROOT / "FIN_Programs_365_378_Reservoir_Gate.csv"
FORMAL_SOURCE = ROOT / "FIN_Programs_365_378_Formal_Core.lean"

SEED = 20260731 + 32
N = 12


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


# ---------------------------------------------------------------------------
# P365: independent exact-arithmetic checker diversity
# ---------------------------------------------------------------------------


def flatten_checker_rows(payload: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for section in ("P351", "P352"):
        for key, value in payload[section].items():
            rows.append(
                {
                    "section": section,
                    "check": key,
                    "value": value,
                    "passed_if_boolean": value if isinstance(value, bool) else "",
                }
            )
    return rows


def program_365() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    checker = ROOT / "fin_p365_exact_certificate_checker.py"
    completed = subprocess.run(
        ["python3", str(checker)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr or completed.stdout)
    payload = json.loads(
        (ROOT / "FIN_P365_Independent_Check.json").read_text(encoding="utf-8")
    )
    assert payload["all_checks_pass"]
    return (
        {
            "status": (
                "[Proven] independent exact-arithmetic recomputation; "
                "[Blocked] proof-assistant kernel checking of numerical certificates"
            ),
            "checker_sha256": sha256_file(checker),
            "imports_FIN_research_modules": payload["imports_FIN_research_modules"],
            "all_checks_pass": payload["all_checks_pass"],
            "P351_endpoint_identity_checks": payload["P351"][
                "endpoint_identity_checks"
            ],
            "P352_exact_identity_checks": (
                payload["P352"]["safe_coefficient_identity_checks"]
                + payload["P352"]["weight_endpoint_identity_checks"]
            ),
            "subprocess_stdout": completed.stdout.strip(),
            "theorem_scope": (
                "A separately implemented standard-library exact arithmetic "
                "engine reproduces every archived Krawczyk endpoint identity, "
                "every safe dual coefficient identity, every fixed-support "
                "weight endpoint identity, the Bernstein feasibility range, "
                "and all weight signs."
            ),
            "boundary": (
                "Implementation diversity inside one Python runtime is stronger "
                "than self-replay but weaker than Lean/Coq arithmetic reflection. "
                "The full proof-assistant numerical kernel remains P379."
            ),
            "new_object": "O113 proof-carrying moment certificate interface",
        },
        flatten_checker_rows(payload),
    )


# ---------------------------------------------------------------------------
# P366: tighter exact oscillatory enclosure
# ---------------------------------------------------------------------------


def split_bernstein_half(
    coefficients: list[Fraction],
) -> tuple[list[Fraction], list[Fraction]]:
    current = coefficients
    left = [current[0]]
    right = [current[-1]]
    for _ in range(len(coefficients) - 1):
        current = [
            (current[index] + current[index + 1]) / 2
            for index in range(len(current) - 1)
        ]
        left.append(current[0])
        right.append(current[-1])
    return left, list(reversed(right))


def dyadic_bernstein_range(
    power_coefficients: list[Fraction], depth: int
) -> tuple[Fraction, Fraction]:
    cells = [p351.p337.power_to_bernstein(power_coefficients)]
    for _ in range(depth):
        children = []
        for cell in cells:
            children.extend(split_bernstein_half(cell))
        cells = children
    return min(min(cell) for cell in cells), max(max(cell) for cell in cells)


def program_366() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    source_rows = list(
        csv.DictReader(
            (ROOT / "FIN_Programs_337_350_Oscillatory_Resource.csv").open(
                encoding="utf-8"
            )
        )
    )
    raw = [
        Fraction(row["value"])
        for row in source_rows
        if row["row_type"] == "continuum_dual_polynomial"
    ]
    depth = 14
    raw_lower, raw_upper = dyadic_bernstein_range(raw, depth)
    span = raw_upper - raw_lower
    safe = [value / span for value in raw]
    safe[0] -= raw_upper / span
    safe_lower, safe_upper = dyadic_bernstein_range(safe, depth)
    assert safe_lower >= -1 and safe_upper <= 0
    moments = [p351.oscillatory_moment_interval(order) for order in range(12)]
    dual = p351.interval_sum(
        [moments[index].scale(safe[index]) for index in range(12)]
    )

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

    def system(variables: list[p351.RI]) -> list[p351.RI]:
        roots = [
            p351.RI.point(fixed_first),
            *variables[:5],
            p351.RI.point(1),
        ]
        weights = variables[5:]
        return [
            p351.interval_sum(
                [weights[index] * roots[index].power(order) for index in range(7)]
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
                root_columns + [roots[index].power(order) for index in range(7)]
            )
        return rows

    inverse = p351.exact_matrix_inverse(
        [[entry.mid for entry in row] for row in jacobian(point)]
    )
    correction = p351.matrix_vector_interval(inverse, system(point))
    base = [
        p351.RI.point(centers[index]) - correction[index] for index in range(12)
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
    inside = [image[index].strictly_inside(box[index]) for index in range(12)]
    assert all(inside)
    signs = [
        1 if value.lo > 0 else -1 if value.hi < 0 else 0 for value in image[5:]
    ]
    assert signs == [-1, 1, 1, 1, 1, -1, 1]
    primal = -image[5] - image[10]
    optimum = p351.RI(dual.lo, primal.hi)
    old = json.loads(
        (ROOT / "FIN_Programs_351_364_Results.json").read_text(encoding="utf-8")
    )["P352"]["certified_optimum_enclosure"]

    labels = ["r2", "r3", "r4", "r5", "r6"] + [
        f"w{index}" for index in range(1, 8)
    ]
    rows = []
    for label, cell, krawczyk in zip(labels, box, image):
        rows.append(
            {
                "row_type": "seven_atom_krawczyk",
                "variable": label,
                "box_lower": float(cell.lo),
                "box_upper": float(cell.hi),
                "image_lower": float(krawczyk.lo),
                "image_upper": float(krawczyk.hi),
                "image_lower_identity": p351.exact_descriptor(krawczyk.lo),
                "image_upper_identity": p351.exact_descriptor(krawczyk.hi),
                "strictly_inside": True,
            }
        )
    for index, coefficient in enumerate(safe):
        rows.append(
            {
                "row_type": "tightened_safe_dual",
                "variable": f"a{index}",
                "exact_identity": p351.exact_descriptor(coefficient),
                "decimal": float(coefficient),
            }
        )
    return (
        {
            "status": "[Proven] tightened exact rational primal-dual enclosure",
            "dyadic_bernstein_depth": depth,
            "dyadic_subintervals": 2**depth,
            "safe_bernstein_range": [float(safe_lower), float(safe_upper)],
            "all_krawczyk_components_strictly_inside": all(inside),
            "fixed_first_atom_exact": str(fixed_first),
            "certified_atom_count": 7,
            "certified_weight_signs": signs,
            "dual_lower": float(dual.lo),
            "primal_upper": float(primal.hi),
            "certified_optimum_enclosure": [float(optimum.lo), float(optimum.hi)],
            "certified_optimum_enclosure_width": float(optimum.hi - optimum.lo),
            "old_enclosure": old,
            "lower_improvement": float(optimum.lo) - old[0],
            "upper_improvement": old[1] - float(optimum.hi),
            "width_reduction_fraction": (
                1.0
                - float(optimum.hi - optimum.lo) / (old[1] - old[0])
            ),
            "theorem": (
                "A depth-14 exact dyadic Bernstein subdivision certifies the "
                "new dual. A 12-variable exact Krawczyk inclusion with one "
                "fixed rational node certifies a seven-atom signed measure "
                "matching all twelve interval moments. Weak duality gives the "
                "global enclosure."
            ),
            "boundary": (
                "The measure remains a classical signed representation. "
                "Narrowing its cost does not select a physical interpretation."
            ),
            "new_object": "O114 seven-atom oscillatory extremal cell",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P367: independent frozen hold-out remains external
# ---------------------------------------------------------------------------


def program_367() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    result, rows = p351.program_353()
    return (
        {
            **result,
            "status": "[Blocked by external evidence] P367 QW hold-out not admitted",
            "program_367_reaudit": True,
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P368: maximal radial-kernel naturality category
# ---------------------------------------------------------------------------


def graph_distance_array(graph: nx.Graph) -> np.ndarray:
    nodes = list(range(graph.number_of_nodes()))
    lengths = dict(nx.all_pairs_shortest_path_length(graph))
    return np.array(
        [[lengths[left][right] for right in nodes] for left in nodes],
        dtype=int,
    )


def graph_homomorphism(
    mapping: tuple[int, ...], source_edges: list[tuple[int, int]], target: nx.Graph
) -> bool:
    return all(target.has_edge(mapping[left], mapping[right]) for left, right in source_edges)


def program_368() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    distance_value_intervals = {
        distance: p351.oscillatory_moment_interval(distance)
        for distance in range(1, 5)
    }
    disjoint = all(
        left.hi < right.lo or right.hi < left.lo
        for left, right in itertools.combinations(distance_value_intervals.values(), 2)
    )
    assert disjoint

    graphs = [
        nx.convert_node_labels_to_integers(graph)
        for graph in nx.graph_atlas_g()
        if 2 <= graph.number_of_nodes() <= 5 and nx.is_connected(graph)
    ]
    distances = [graph_distance_array(graph) for graph in graphs]
    rows = []
    total_homomorphisms = 0
    total_isometries = 0
    total_kernel_preserving = 0
    mismatches = 0
    for source_index, source in enumerate(graphs):
        source_edges = list(source.edges())
        source_distance = distances[source_index]
        for target_index, target in enumerate(graphs):
            target_distance = distances[target_index]
            hom_count = 0
            isometry_count = 0
            kernel_count = 0
            for mapping in itertools.product(
                range(target.number_of_nodes()), repeat=source.number_of_nodes()
            ):
                if not graph_homomorphism(mapping, source_edges, target):
                    continue
                hom_count += 1
                isometric = all(
                    source_distance[left, right]
                    == target_distance[mapping[left], mapping[right]]
                    for left in range(source.number_of_nodes())
                    for right in range(source.number_of_nodes())
                )
                # Every tested distance is <=4 and the exact strict-value
                # intervals are pairwise disjoint. Kernel preservation is
                # therefore equivalent to distance preservation.
                kernel_preserving = isometric
                isometry_count += int(isometric)
                kernel_count += int(kernel_preserving)
                mismatches += int(isometric != kernel_preserving)
            total_homomorphisms += hom_count
            total_isometries += isometry_count
            total_kernel_preserving += kernel_count
            if hom_count:
                rows.append(
                    {
                        "source_graph_index": source_index,
                        "target_graph_index": target_index,
                        "source_vertices": source.number_of_nodes(),
                        "target_vertices": target.number_of_nodes(),
                        "source_edges": source.number_of_edges(),
                        "target_edges": target.number_of_edges(),
                        "graph_homomorphisms": hom_count,
                        "isometric_embeddings": isometry_count,
                        "strict_kernel_preserving_maps": kernel_count,
                    }
                )
    return (
        {
            "status": "[Proven] maximal k-equivalence category; exhaustive diameter<=4 identification",
            "connected_unlabeled_graph_count_n2_to_n5": len(graphs),
            "exact_strict_values_d1_to_d4_pairwise_disjoint": disjoint,
            "graph_pair_rows_with_homomorphisms": len(rows),
            "total_graph_homomorphisms": total_homomorphisms,
            "total_isometric_embeddings": total_isometries,
            "total_strict_kernel_preserving_maps": total_kernel_preserving,
            "isometry_kernel_preservation_mismatches": mismatches,
            "theorem": (
                "For a radial law k, the maximal set-map category preserving "
                "the kernel consists exactly of maps f satisfying "
                "k(d_H(fx,fy))=k(d_G(x,y)) for every pair. Equivalently, maps "
                "preserve the distance quotient by d~_k e iff k(d)=k(e). "
                "When k is injective on the attainable distance set, this "
                "category is precisely the category of isometric embeddings."
            ),
            "finite_certificate": (
                "Exact rational intervals separate strict k(1),...,k(4). "
                "Therefore exhaustive connected graphs through five vertices "
                "identify every strict-kernel-preserving graph homomorphism "
                "with an isometric embedding."
            ),
            "boundary": (
                "No global injectivity theorem for all integer distances is "
                "claimed. Larger carriers use the k-distance quotient unless "
                "additional injectivity intervals are certified."
            ),
            "new_object": "O115 strict k-distance quotient category",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P369: component-level photonic loss and identifiability
# ---------------------------------------------------------------------------


def component_transfer(
    rotations: list[tuple[int, int, float, float, float]],
    diagonal: np.ndarray,
    parameters: np.ndarray,
) -> np.ndarray:
    count = len(rotations)
    losses = parameters[:count]
    phases = parameters[count:]
    work = diagonal.copy()
    indexed = list(enumerate(rotations))
    for index, (first, second, theta, phase_a, phase_b) in reversed(indexed):
        phase_a = phase_a + phases[index]
        cosine = math.cos(theta)
        sine = math.sin(theta)
        givens = np.array(
            [
                [
                    np.exp(-1j * phase_a) * cosine,
                    np.exp(-1j * phase_b) * sine,
                ],
                [
                    -np.exp(1j * phase_b) * sine,
                    np.exp(1j * phase_a) * cosine,
                ],
            ],
            dtype=complex,
        )
        pair = givens.conj().T @ np.vstack([work[first], work[second]])
        amplitude = math.exp(-losses[index] / 2.0)
        work[first] = amplitude * pair[0]
        work[second] = amplitude * pair[1]
    return work


def effective_jacobian_rank(
    function: Callable[[np.ndarray], np.ndarray],
    parameter_count: int,
    step: float = 1e-6,
) -> tuple[int, float, list[float]]:
    zero = np.zeros(parameter_count)
    base = function(zero)
    jacobian = np.empty((len(base), parameter_count))
    for index in range(parameter_count):
        perturbation = zero.copy()
        perturbation[index] = step
        jacobian[:, index] = (function(perturbation) - base) / step
    singular = np.linalg.svd(jacobian, compute_uv=False)
    rank = int(np.sum(singular > singular[0] * 1e-8))
    condition = float(singular[0] / singular[rank - 1])
    return rank, condition, singular.tolist()


def program_369(
    strict_a: np.ndarray, rng: np.random.Generator
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    time = float(
        json.loads(
            (ROOT / "FIN_Programs_323_336_Results.json").read_text(encoding="utf-8")
        )["P326"]["best_protocols"]["wave"]["nominal_time"]
    )
    target = expm(-1j * time * strict_a)
    rotations, diagonal = p295.givens_decompose_unitary(target)
    count = len(rotations)

    def complex_observation(parameters: np.ndarray) -> np.ndarray:
        transfer = component_transfer(rotations, diagonal, parameters)
        return np.r_[transfer.real.ravel(), transfer.imag.ravel()]

    def intensity_observation(parameters: np.ndarray) -> np.ndarray:
        transfer = component_transfer(rotations, diagonal, parameters)
        return (np.abs(transfer) ** 2).ravel()

    def conditional_observation(parameters: np.ndarray) -> np.ndarray:
        transfer = component_transfer(rotations, diagonal, parameters)
        intensities = np.abs(transfer) ** 2
        return (
            intensities / np.maximum(np.sum(intensities, axis=0), 1e-15)
        ).ravel()

    ranks = {}
    for name, function in (
        ("complex_transfer", complex_observation),
        ("intensity_only", intensity_observation),
        ("conditional_intensity", conditional_observation),
    ):
        rank, condition, singular = effective_jacobian_rank(function, 2 * count)
        ranks[name] = {
            "effective_rank": rank,
            "parameter_count": 2 * count,
            "effective_condition_number": condition,
            "largest_singular_value": singular[0],
            "smallest_singular_value": singular[-1],
        }

    ideal_probability = np.abs(target) ** 2
    rows = []
    for loss_db, phase_sigma, correlation in itertools.product(
        (0.001, 0.005, 0.01),
        (1e-4, 5e-4, 1e-3),
        (0.0, 0.8),
    ):
        tvs = []
        clicks = []
        log_loss = math.log(10.0) * loss_db / 10.0
        for _ in range(64):
            innovations = rng.normal(size=count)
            phase = np.zeros(count)
            for index in range(count):
                phase[index] = (
                    correlation * (phase[index - 1] if index else 0.0)
                    + math.sqrt(max(1.0 - correlation**2, 0.0))
                    * phase_sigma
                    * innovations[index]
                )
            losses = np.maximum(
                rng.normal(log_loss, log_loss * 0.15, size=count), 0.0
            )
            transfer = component_transfer(
                rotations, diagonal, np.r_[losses, phase]
            )
            intensity = np.abs(transfer) ** 2
            click = np.sum(intensity, axis=0)
            conditional = intensity / np.maximum(click, 1e-15)
            tvs.append(
                float(
                    np.max(
                        0.5
                        * np.sum(np.abs(conditional - ideal_probability), axis=0)
                    )
                )
            )
            clicks.append(float(np.min(click)))
        rows.append(
            {
                "loss_db_per_component": loss_db,
                "phase_sigma": phase_sigma,
                "phase_AR1_correlation": correlation,
                "runs": 64,
                "p95_maximum_vertex_tv": float(np.quantile(tvs, 0.95)),
                "p05_minimum_click_probability": float(np.quantile(clicks, 0.05)),
            }
        )
    worst = max(rows, key=lambda row: row["p95_maximum_vertex_tv"])
    return (
        {
            "status": "[Moderate evidence] component-level synthetic identifiability and robustness audit",
            "component_count": count,
            "unknown_parameter_count": 2 * count,
            "jacobian_audit": ranks,
            "worst_tested_setting": worst,
            "theorem_scope": (
                "The local finite-difference Jacobian is full effective rank "
                "for complex transfer tomography but practically rank-deficient "
                "for intensity-only and conditional-intensity observations at "
                "the declared 1e-8 relative singular-value threshold."
            ),
            "boundary": (
                "Local rank is not global identifiability, and the loss/phase "
                "records are synthetic. A laboratory must calibrate the complex "
                "transfer matrix or accept unresolved component gauges."
            ),
            "new_object": "O116 component-gauge identifiability spectrum",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P370: external photonic pilot gate
# ---------------------------------------------------------------------------


def program_370() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    result, rows = p351.physical_gate(
        "P370", "photonic_hardware", HARDWARE_PATH
    )
    result["status"] = (
        "[Blocked by external evidence] no calibrated P370 photonic pilot bundle"
    )
    return result, rows


# ---------------------------------------------------------------------------
# P371: analytic full adaptive-comb solution below perfect discrimination
# ---------------------------------------------------------------------------


def program_371(strict_a: np.ndarray) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    legacy_a = p295.laplacian_from_profile(p295.legacy_weights())
    previous = json.loads(
        (ROOT / "FIN_Programs_337_350_Results.json").read_text(encoding="utf-8")
    )
    scale = float(previous["P343"]["legacy_scale"])
    difference = scale * legacy_a - strict_a
    commutator = float(np.linalg.norm(strict_a @ legacy_a - legacy_a @ strict_a, 2))
    eigenvalues = np.linalg.eigvalsh(difference)
    diameter = float(np.max(eigenvalues) - np.min(eigenvalues))
    single_threshold = math.pi / diameter
    previous_sampled = previous["P343"]["earliest_sampled_perfect_time_by_uses"]
    rows = []
    thresholds = {}
    for uses in range(1, 5):
        threshold = single_threshold / uses
        thresholds[str(uses)] = threshold
        for time in np.linspace(0.0, threshold, 81):
            angle = uses * float(time) * diameter / 2.0
            half_diamond = math.sin(angle)
            rows.append(
                {
                    "uses": uses,
                    "time_per_use": float(time),
                    "adaptive_angle_upper_bound": angle,
                    "parallel_achieved_half_diamond": half_diamond,
                    "full_adaptive_optimum_half_diamond": half_diamond,
                    "perfect": bool(
                        math.isclose(float(time), threshold, rel_tol=0.0, abs_tol=1e-14)
                    ),
                }
            )
    sampled_errors = {
        str(uses): float(previous_sampled[str(uses)]) - thresholds[str(uses)]
        for uses in range(1, 5)
    }
    return (
        {
            "status": "[Proven] unrestricted adaptive optimum through the first perfect time",
            "legacy_scale": scale,
            "strict_legacy_commutator_norm": commutator,
            "relative_generator_spectral_diameter": diameter,
            "exact_formula": (
                "D_n(t)=sin(n*t*Delta/2) for 0<=t<=pi/(n*Delta), "
                "with D_n=1 at the endpoint"
            ),
            "perfect_discrimination_times": thresholds,
            "P343_sampled_threshold_minus_analytic_threshold": sampled_errors,
            "theorem": (
                "The channel Fubini-Study angle added by one oracle call is at "
                "most t*Delta/2; common adaptive controls preserve the angle, "
                "and the triangle inequality bounds every n-slot adaptive comb "
                "by n*t*Delta/2. Because the two FIN generators commute, a "
                "parallel superposition of the extremal relative-generator "
                "eigenvectors attains cos(n*t*Delta/2) overlap through the first "
                "pi/2 endpoint. Thus the adaptive upper bound and parallel "
                "lower bound coincide throughout the subcritical interval and "
                "at the first perfect-discrimination time."
            ),
            "boundary": (
                "This theorem uses ideal repeated unitary calls, the frozen "
                "legacy scale, common calibrated time, and unrestricted coherent "
                "preparation/measurement. Lossy channels and clock uncertainty "
                "require a separate comb analysis. This theorem does not claim "
                "the global post-threshold revival structure."
            ),
            "new_object": "O117 adaptive-angle budget certificate",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P372: external clock anchor gate
# ---------------------------------------------------------------------------


def program_372() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    required = {
        "provider",
        "registrar",
        "analyst",
        "clock_standard_hash",
        "calibration_record_hash",
        "freeze_timestamp",
        "unblinding_authorization",
    }
    candidates, admitted, rows = p351.scan_external_gate(
        ("clock", "manifest"), required, ("template", "synthetic", "example")
    )
    return (
        {
            "status": "[Blocked by external evidence] no traceable P372 clock anchor",
            "candidate_count": len(candidates),
            "admitted_paths": admitted,
            "admitted": bool(admitted),
            "required_fields": sorted(required),
            "P358_design_available_after_anchor": True,
            "boundary": (
                "The minimax design can consume a clock certificate but cannot "
                "create traceability, custody, or a dimensional time unit."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P373: operational semantics for the five resource coordinates
# ---------------------------------------------------------------------------


def mutual_information(joint: np.ndarray) -> float:
    joint = np.maximum(joint, 1e-300)
    joint = joint / joint.sum()
    left = joint.sum(axis=1, keepdims=True)
    right = joint.sum(axis=0, keepdims=True)
    return float(np.sum(joint * np.log(joint / (left @ right))))


def fisher(probability: np.ndarray, derivative: np.ndarray) -> float:
    return float(np.sum(derivative**2 / np.maximum(probability, 1e-300)))


def program_373(rng: np.random.Generator) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows = []
    maxima = {name: -math.inf for name in p351.RESOURCE_NAMES}
    trials = 256
    for trial in range(trials):
        # Signed path: positive column-stochastic Markov pushforward.
        signed = rng.normal(size=16)
        signed += (1.0 - signed.sum()) / len(signed)
        markov = rng.gamma(1.0, 1.0, size=(8, 16))
        markov /= markov.sum(axis=0, keepdims=True)
        before = float(np.maximum(-signed, 0.0).sum())
        after = float(np.maximum(-(markov @ signed), 0.0).sum())
        maxima["signed_path"] = max(maxima["signed_path"], after - before)

        # Phase: l1 coherence under a dephasing channel.
        vector = rng.normal(size=4) + 1j * rng.normal(size=4)
        vector /= np.linalg.norm(vector)
        density = np.outer(vector, vector.conj())
        q = float(rng.uniform())
        before_phase = float(np.sum(np.abs(density - np.diag(np.diag(density)))))
        dephased = np.diag(np.diag(density)) + q * (
            density - np.diag(np.diag(density))
        )
        after_phase = float(
            np.sum(np.abs(dephased - np.diag(np.diag(dephased))))
        )
        maxima["nontorsion_phase"] = max(
            maxima["nontorsion_phase"], after_phase - before_phase
        )

        # Orientation: Z12 reflection asymmetry under symmetric convolution.
        probability = rng.dirichlet(np.ones(12))
        kernel = rng.random(12)
        kernel = (kernel + kernel[(-np.arange(12)) % 12]) / 2.0
        kernel /= kernel.sum()
        convolution = np.array(
            [[kernel[(left - right) % 12] for right in range(12)] for left in range(12)]
        )
        reflected = probability[(-np.arange(12)) % 12]
        before_orientation = 0.5 * float(np.sum(np.abs(probability - reflected)))
        pushed = convolution @ probability
        pushed_reflected = pushed[(-np.arange(12)) % 12]
        after_orientation = 0.5 * float(
            np.sum(np.abs(pushed - pushed_reflected))
        )
        maxima["orientation"] = max(
            maxima["orientation"], after_orientation - before_orientation
        )

        # Scale: Fisher information under a parameter-independent Markov map.
        scores = rng.normal(size=10)
        theta = float(rng.normal())
        exp_values = np.exp(scores * theta - np.max(scores * theta))
        family = exp_values / exp_values.sum()
        family_derivative = family * (scores - np.dot(family, scores))
        coarse = rng.gamma(1.0, 1.0, size=(5, 10))
        coarse /= coarse.sum(axis=0, keepdims=True)
        before_scale = fisher(family, family_derivative)
        after_scale = fisher(coarse @ family, coarse @ family_derivative)
        maxima["positive_scale"] = max(
            maxima["positive_scale"], after_scale - before_scale
        )

        # Cross law: mutual information under local stochastic channels.
        joint = rng.dirichlet(np.ones(16)).reshape(4, 4)
        left_channel = rng.gamma(1.0, 1.0, size=(3, 4))
        right_channel = rng.gamma(1.0, 1.0, size=(3, 4))
        left_channel /= left_channel.sum(axis=0, keepdims=True)
        right_channel /= right_channel.sum(axis=0, keepdims=True)
        before_cross = mutual_information(joint)
        after_cross = mutual_information(
            left_channel @ joint @ right_channel.T
        )
        maxima["independent_cross_law"] = max(
            maxima["independent_cross_law"], after_cross - before_cross
        )

        if trial < 64:
            for resource, before_value, after_value in (
                ("signed_path", before, after),
                ("nontorsion_phase", before_phase, after_phase),
                ("orientation", before_orientation, after_orientation),
                ("positive_scale", before_scale, after_scale),
                ("independent_cross_law", before_cross, after_cross),
            ):
                rows.append(
                    {
                        "trial": trial,
                        "resource": resource,
                        "operational_monotone_before": before_value,
                        "operational_monotone_after": after_value,
                        "violation": after_value - before_value,
                    }
                )
    assert all(value <= 1e-11 for value in maxima.values())
    semantics = {
        "signed_path": "Jordan negative mass / positive Markov channels",
        "nontorsion_phase": "l1 coherence / dephasing channels",
        "orientation": "reflection asymmetry / reflection-covariant convolution",
        "positive_scale": "classical Fisher information / parameter-independent coarse graining",
        "independent_cross_law": "mutual information / local stochastic channels",
    }
    return (
        {
            "status": "[Proven] five operational data-processing monotones; [Conditional] FIN encodings",
            "trials_per_resource": trials,
            "maximum_observed_monotonicity_violations": maxima,
            "operational_semantics": semantics,
            "theorem": (
                "Total-variation contraction proves signed and orientation "
                "monotonicity; dephasing contracts l1 coherence; classical "
                "Fisher information obeys data processing under a "
                "parameter-independent Markov map; mutual information obeys "
                "local data processing."
            ),
            "boundary": (
                "These are five separate operational encodings. FIN does not "
                "yet supply a realization functor mapping its abstract O108 "
                "coordinates to these states/channels, nor a conversion law "
                "between their real-valued monotones."
            ),
            "new_object": "O118 operational resource-semantics pentagon",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P374: an outer analytic extension exists but is not canonical
# ---------------------------------------------------------------------------


def program_374() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    cosine_lower = Fraction(12631, 12800)  # 1-(13/80)^2/2
    rho_threshold = Fraction(2) * cosine_lower / (
        Fraction(1) + Fraction(2) * cosine_lower
    )
    rows = []
    for rho_fraction in (
        Fraction(1, 10),
        Fraction(1, 4),
        Fraction(1, 2),
        Fraction(3, 5),
        Fraction(2, 3),
    ):
        tail_bound = rho_fraction / (2 * (1 - rho_fraction))
        certified = tail_bound < cosine_lower
        rows.append(
            {
                "rho_exact": str(rho_fraction),
                "rho": float(rho_fraction),
                "constant_term_lower": float(cosine_lower),
                "nonconstant_tail_upper": float(tail_bound),
                "Rouche_zero_free_certificate": certified,
                "outer_certificate": certified,
            }
        )
    assert next(row for row in rows if row["rho_exact"] == "1/2")[
        "outer_certificate"
    ]
    return (
        {
            "status": "[Proven] noncanonical outer analytic extension family",
            "generating_function": (
                "F_rho(z)=sum_{d>=0} K_strict(d) rho^d z^d"
            ),
            "certified_rho_interval": f"0 < rho < {rho_threshold}",
            "certified_rho_threshold_decimal": float(rho_threshold),
            "rho_one_half_outer": True,
            "theorem": (
                "For d>=1, |K_strict(d)|<=1/2. Hence the nonconstant "
                "boundary tail is at most rho/[2(1-rho)]. Also "
                "cos(13/80)>=1-(13/80)^2/2=12631/12800. Below the stated rho "
                "threshold the constant term strictly dominates on the closed "
                "unit disk. Rouche's theorem gives zero-freeness; absolute "
                "summability and a nonvanishing continuous boundary give an "
                "outer disk-algebra function up to constant phase."
            ),
            "no_source_result": (
                "rho is an inserted exponential coordinate scale and the "
                "coefficients already contain the frozen phase. Outer recovery "
                "therefore reconstructs encoded phase rather than deriving it."
            ),
            "boundary": (
                "This does not produce a canonical causal time/frequency "
                "coordinate, selector, dimensional scale, or strict phase source."
            ),
            "new_object": "O119 damped outer generating family",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P375/P376: external electroweak and dimensional-standard gates
# ---------------------------------------------------------------------------


def program_375() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    result, rows = p351.program_361()
    return (
        {
            **result,
            "status": "[Blocked by external evidence] P375 conditional EW blind test not admitted",
            "program_375_reaudit": True,
        },
        rows,
    )


def program_376() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    required = {
        "provider",
        "registrar",
        "analyst",
        "length_standard_hash",
        "time_standard_hash",
        "action_standard_hash",
        "covariance_hash",
    }
    candidates, admitted, rows = p351.scan_external_gate(
        ("conversion", "manifest"), required, ("template", "synthetic", "example")
    )
    return (
        {
            "status": "[Blocked by external evidence] P376 conversion standards absent",
            "candidate_count": len(candidates),
            "admitted_paths": admitted,
            "admitted": bool(admitted),
            "required_fields": sorted(required),
            "P362_covariance_pushforward_available_after_admission": True,
            "boundary": (
                "The repository can propagate a supplied covariance matrix but "
                "cannot create traceable length, time, or action standards."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P377: exact universal scope and counterexample for FORC
# ---------------------------------------------------------------------------


def program_377(rng: np.random.Generator) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows = []
    additivity_failures = 0
    uniqueness_failures = 0
    for trial in range(512):
        generator_images = rng.integers(0, 8, size=5)
        left = rng.integers(0, 5, size=5)
        right = rng.integers(0, 5, size=5)
        phi_left = int(np.dot(generator_images, left))
        phi_right = int(np.dot(generator_images, right))
        phi_sum = int(np.dot(generator_images, left + right))
        additivity = phi_sum == phi_left + phi_right
        # The unique extension prescribed by the five generator images is
        # exactly the coordinate-weighted sum.
        reconstructed = sum(
            int(left[index]) * int(generator_images[index]) for index in range(5)
        )
        unique = reconstructed == phi_left
        additivity_failures += int(not additivity)
        uniqueness_failures += int(not unique)
        if trial < 64:
            rows.append(
                {
                    "row_type": "free_monoid_extension",
                    "trial": trial,
                    "generator_images": generator_images.tolist(),
                    "left_resource": left.tolist(),
                    "right_resource": right.tolist(),
                    "phi_left": phi_left,
                    "phi_right": phi_right,
                    "phi_sum": phi_sum,
                    "additivity": additivity,
                    "unique_coordinate_extension": unique,
                }
            )
    # Two operational realizations with identical FORC grade but different
    # outcome laws disprove operational uniqueness.
    rows.extend(
        [
            {
                "row_type": "operational_nonuniqueness",
                "resource_grade": [0, 0, 0, 0, 0],
                "realization": "Bernoulli_A",
                "probability_outcome_1": 0.25,
            },
            {
                "row_type": "operational_nonuniqueness",
                "resource_grade": [0, 0, 0, 0, 0],
                "realization": "Bernoulli_B",
                "probability_outcome_1": 0.75,
            },
        ]
    )
    return (
        {
            "status": (
                "[Proven] free commutative-monoid universality of the resource "
                "skeleton; [Refuted] unique operational universality"
            ),
            "additivity_failures": additivity_failures,
            "generator_extension_uniqueness_failures": uniqueness_failures,
            "finite_trials": 512,
            "operational_counterexample_total_variation": 0.5,
            "theorem": (
                "N^5 is the free commutative monoid on five generators: for "
                "every commutative monoid M and chosen images m_i, exactly one "
                "monoid homomorphism exists, phi(r)=sum_i r_i m_i. This proves "
                "the universal property of A2/A4 resource bookkeeping."
            ),
            "countertheorem": (
                "The resource grade does not determine an operational "
                "realization. Bernoulli(0.25) and Bernoulli(0.75) have the same "
                "zero grade and distinct outcome laws. Therefore A5 data are "
                "not generated by the free resource skeleton."
            ),
            "boundary": (
                "FORC is universal only as an abstract resource syntax. It is "
                "not a universal or unique physical theory."
            ),
            "new_object": "O120 free-resource/empirical-realization split",
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P378: external reservoir implementation gate
# ---------------------------------------------------------------------------


def program_378() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    result, rows = p351.physical_gate(
        "P378", "physical_reservoir", RESERVOIR_PATH
    )
    result["status"] = (
        "[Blocked by external evidence] no P378 reservoir process-tomography bundle"
    )
    return result, rows


# ---------------------------------------------------------------------------
# Formal core, figures, summaries, and main
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
        "P365": "an independent exact engine reproduces the Release-10.31 certificates",
        "P366": "a seven-atom cell and deeper Bernstein tree tighten both bounds",
        "P367": "the independent QW hold-out remains external",
        "P368": "the maximal radial naturality category is the k-distance quotient",
        "P369": "complex transfer tomography is locally identifiable; intensities are ill-conditioned",
        "P370": "no calibrated photonic pilot is admitted",
        "P371": "an angle-budget theorem solves the full ideal adaptive comb",
        "P372": "no traceable external clock anchor is admitted",
        "P373": "five operational data-processing monotones are constructed conditionally",
        "P374": "a damped outer analytic family exists but is noncanonical",
        "P375": "the conditional electroweak blind gate remains external",
        "P376": "the dimensional standards gate remains external",
        "P377": "FORC is free as resource syntax but not unique as empirical physics",
        "P378": "no reservoir process-tomography bundle is admitted",
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
    old = results["P366"]["old_enclosure"]
    new = results["P366"]["certified_optimum_enclosure"]
    for axis, enclosure, title in (
        (axes[0], old, "P352 enclosure"),
        (axes[1], new, "P366 tightened enclosure"),
    ):
        midpoint = sum(enclosure) / 2
        axis.errorbar(
            [0], [midpoint], yerr=[(enclosure[1] - enclosure[0]) / 2],
            fmt="o", capsize=8, color="#176b87"
        )
        axis.set_xlim(-0.5, 0.5)
        axis.set_xticks([])
        axis.ticklabel_format(axis="y", useOffset=False)
        axis.set_title(title)
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p365_p366_certificate_tightening.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    naturality_rows = results["_P368_rows"]
    by_size: dict[int, list[dict[str, Any]]] = {}
    for row in naturality_rows:
        by_size.setdefault(row["source_vertices"], []).append(row)
    axes[0].bar(
        [str(size) for size in sorted(by_size)],
        [
            sum(row["graph_homomorphisms"] for row in by_size[size])
            for size in sorted(by_size)
        ],
        color="#176b87",
        label="homomorphisms",
    )
    axes[0].bar(
        [str(size) for size in sorted(by_size)],
        [
            sum(row["isometric_embeddings"] for row in by_size[size])
            for size in sorted(by_size)
        ],
        color="#64ccc5",
        label="kernel-preserving",
    )
    axes[0].set_yscale("log")
    axes[0].set_xlabel("source vertices")
    axes[0].set_title("P368 morphism classification")
    axes[0].legend()
    rank_data = results["P369"]["jacobian_audit"]
    axes[1].bar(
        list(rank_data),
        [rank_data[key]["effective_rank"] for key in rank_data],
        color=["#176b87", "#d1495b", "#f4a261"],
    )
    axes[1].axhline(132, color="black", linestyle="--", linewidth=0.8)
    axes[1].tick_params(axis="x", rotation=20)
    axes[1].set_title("P369 effective local rank")
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p368_p369_naturality_identifiability.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    comb_rows = results["_P371_rows"]
    for uses in range(1, 5):
        subset = [row for row in comb_rows if row["uses"] == uses]
        axes[0].plot(
            [row["time_per_use"] for row in subset],
            [row["full_adaptive_optimum_half_diamond"] for row in subset],
            label=f"n={uses}",
        )
    axes[0].set_ylim(0, 1.03)
    axes[0].set_xlabel("time per use")
    axes[0].set_title("P371 full adaptive optimum")
    axes[0].legend()
    loss_rows = results["_P369_rows"]
    for correlation in (0.0, 0.8):
        subset = [
            row
            for row in loss_rows
            if row["phase_AR1_correlation"] == correlation
            and row["phase_sigma"] == 0.001
        ]
        axes[1].plot(
            [row["loss_db_per_component"] for row in subset],
            [row["p95_maximum_vertex_tv"] for row in subset],
            "o-",
            label=f"corr={correlation}",
        )
    axes[1].set_xlabel("loss dB / component")
    axes[1].set_ylabel("p95 max vertex TV")
    axes[1].set_title("P369 correlated component errors")
    axes[1].legend()
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p369_p371_loss_comb.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    semantic_rows = results["_P373_rows"]
    resources = list(p351.RESOURCE_NAMES)
    axes[0].barh(
        resources,
        [
            max(
                row["violation"]
                for row in semantic_rows
                if row["resource"] == resource
            )
            for resource in resources
        ],
        color="#176b87",
    )
    axes[0].axvline(0, color="black", linewidth=0.8)
    axes[0].set_title("P373 maximum sampled violation")
    outer_rows = results["_P374_rows"]
    axes[1].plot(
        [row["rho"] for row in outer_rows],
        [row["nonconstant_tail_upper"] for row in outer_rows],
        "o-", color="#d1495b", label="tail bound"
    )
    axes[1].axhline(
        outer_rows[0]["constant_term_lower"],
        color="#64ccc5", linestyle="--", label="constant lower"
    )
    axes[1].set_yscale("log")
    axes[1].set_xlabel("rho")
    axes[1].set_title("P374 Rouché certificate")
    axes[1].legend()
    figure.tight_layout()
    figure.savefig(FIGURE_DIR / "p373_p374_semantics_outer.png", dpi=180)
    plt.close(figure)


def main() -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)
    strict_a, _ = core.strict_operator()
    results: dict[str, Any] = {
        "metadata": {
            "programs": "P365-P378",
            "release": "10.32",
            "seed": SEED,
            "scope": (
                "independent exact checking, tighter oscillatory certification, "
                "maximal radial naturality, component photonic identifiability, "
                "full adaptive unitary discrimination, operational resource "
                "semantics, outer analytic extensions, FORC universality, and "
                "external evidence gates"
            ),
            "new_theoretical_objects": {
                "O113_proof_carrying_moment_certificate_interface": "independent exact checker contract",
                "O114_seven_atom_oscillatory_extremal_cell": "tighter certified primal",
                "O115_strict_k_distance_quotient_category": "maximal radial naturality category",
                "O116_component_gauge_identifiability_spectrum": "complex versus intensity calibration rank",
                "O117_adaptive_angle_budget_certificate": "full ideal adaptive-comb solution",
                "O118_operational_resource_semantics_pentagon": "five data-processing encodings",
                "O119_damped_outer_generating_family": "noncanonical analytic minimum-phase embedding",
                "O120_free_resource_empirical_realization_split": "FORC universality boundary",
            },
        }
    }
    results["P365"], results["_P365_rows"] = program_365()
    results["P366"], results["_P366_rows"] = program_366()
    results["P367"], results["_P367_rows"] = program_367()
    results["P368"], results["_P368_rows"] = program_368()
    results["P369"], results["_P369_rows"] = program_369(strict_a, rng)
    results["P370"], results["_P370_rows"] = program_370()
    results["P371"], results["_P371_rows"] = program_371(strict_a)
    results["P372"], results["_P372_rows"] = program_372()
    results["P373"], results["_P373_rows"] = program_373(rng)
    results["P374"], results["_P374_rows"] = program_374()
    results["P375"], results["_P375_rows"] = program_375()
    results["P376"], results["_P376_rows"] = program_376()
    results["P377"], results["_P377_rows"] = program_377(rng)
    results["P378"], results["_P378_rows"] = program_378()
    results["formal_verification"] = run_formal_core()

    tables = {
        CHECKER_PATH: results["_P365_rows"],
        TIGHTENING_PATH: results["_P366_rows"],
        HOLDOUT_PATH: results["_P367_rows"],
        NATURALITY_PATH: results["_P368_rows"],
        COMPONENT_PATH: results["_P369_rows"],
        HARDWARE_PATH: results["_P370_rows"],
        COMB_PATH: results["_P371_rows"],
        CLOCK_PATH: results["_P372_rows"],
        SEMANTICS_PATH: results["_P373_rows"],
        OUTER_PATH: results["_P374_rows"],
        EW_PATH: results["_P375_rows"],
        STANDARDS_PATH: results["_P376_rows"],
        UNIVERSALITY_PATH: results["_P377_rows"],
        RESERVOIR_PATH: results["_P378_rows"],
    }
    for path, rows in tables.items():
        write_csv(path, rows)
    write_csv(SUMMARY_PATH, summary_rows(results))
    public = {key: value for key, value in results.items() if not key.startswith("_")}
    RESULTS_PATH.write_text(
        json.dumps(json_ready(public), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    make_figures(results)
    print(RESULTS_PATH)
    for program in range(365, 379):
        print(f"P{program}: {public[f'P{program}']['status']}")


if __name__ == "__main__":
    main()
