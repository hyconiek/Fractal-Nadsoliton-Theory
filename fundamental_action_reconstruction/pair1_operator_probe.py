#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
DEFAULT_CONFIG = ROOT / "pair1_operator_probe_config.json"
DEFAULT_REPORT = ROOT / "generated" / "pair1_operator_probe_report.json"


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def resolve_path(relative_path: str) -> Path:
    return REPO / relative_path


def relative_to_repo(path: Path) -> str:
    return str(path.relative_to(REPO))


def get_by_path(payload: Any, path: list[str]) -> Any:
    current = payload
    for key in path:
        if not isinstance(current, dict) or key not in current:
            raise KeyError(".".join(path))
        current = current[key]
    return current


def parse_float(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        text = value.strip()
        try:
            return float(text)
        except ValueError:
            return None
    return None


def matrix_from_rotated_diagonal(lam_parallel: float, lam_perp: float, psi0: float) -> list[list[float]]:
    cpsi = math.cos(psi0)
    spsi = math.sin(psi0)
    a_1 = lam_parallel * cpsi * cpsi + lam_perp * spsi * spsi
    b_1 = (lam_parallel - lam_perp) * cpsi * spsi
    d_1 = lam_parallel * spsi * spsi + lam_perp * cpsi * cpsi
    return [[a_1, b_1], [b_1, d_1]]


def summarize_operator(matrix: list[list[float]], tolerance: float) -> dict[str, Any]:
    a_1 = float(matrix[0][0])
    b_1 = float(matrix[0][1])
    d_1 = float(matrix[1][1])
    delta_1 = [a_1 - d_1, b_1]
    classifier = "TRIVIAL_SELECTOR"
    if abs(delta_1[0]) > tolerance or abs(delta_1[1]) > tolerance:
        classifier = "ANCHOR_IMPORTED_SPLIT"
    return {
        "matrix": matrix,
        "a_1": a_1,
        "b_1": b_1,
        "d_1": d_1,
        "trace_A_1": a_1 + d_1,
        "Delta_1": delta_1,
        "classifier": classifier,
    }


def matrix_difference(lhs: list[list[float]], rhs: list[list[float]]) -> list[list[float]]:
    return [
        [float(lhs[0][0] - rhs[0][0]), float(lhs[0][1] - rhs[0][1])],
        [float(lhs[1][0] - rhs[1][0]), float(lhs[1][1] - rhs[1][1])],
    ]


def load_sources(config: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Path]]:
    cache: dict[str, Any] = {}
    source_paths: dict[str, Path] = {}
    for source_name, relative_path in config["source_paths"].items():
        path = resolve_path(relative_path)
        source_paths[source_name] = path
        if path.suffix == ".json":
            cache[source_name] = load_json(path)
        else:
            cache[source_name] = load_text(path)
    return cache, source_paths


def collect_required_inputs(
    config: dict[str, Any], sources: dict[str, Any], source_paths: dict[str, Path]
) -> tuple[dict[str, dict[str, Any]], list[str]]:
    values: dict[str, dict[str, Any]] = {}
    missing: list[str] = []
    for item in config["required_numeric_inputs"]:
        source_name = item["source"]
        field_path = item["path"]
        path_label = f"{relative_to_repo(source_paths[source_name])}::" + ".".join(field_path)
        try:
            raw_value = get_by_path(sources[source_name], field_path)
        except KeyError:
            missing.append(path_label)
            continue
        numeric = parse_float(raw_value)
        if numeric is None:
            missing.append(path_label)
            continue
        values[item["id"]] = {
            "value": numeric,
            "source": source_name,
            "path": path_label,
        }
    return values, missing


def run_frontier_checks(
    config: dict[str, Any], sources: dict[str, Any], source_paths: dict[str, Path]
) -> tuple[list[dict[str, Any]], list[str]]:
    checks: list[dict[str, Any]] = []
    missing: list[str] = []
    for item in config["frozen_frontier_checks"]:
        source_name = item["source"]
        field_path = item["path"]
        expect = item["equals"]
        path_label = f"{relative_to_repo(source_paths[source_name])}::" + ".".join(field_path)
        try:
            actual = get_by_path(sources[source_name], field_path)
        except KeyError:
            checks.append(
                {
                    "id": item["id"],
                    "status": "MISSING",
                    "path": path_label,
                    "expected": expect,
                }
            )
            missing.append(path_label)
            continue
        ok = actual == expect
        checks.append(
            {
                "id": item["id"],
                "status": "PASS" if ok else "FAIL",
                "path": path_label,
                "actual": actual,
                "expected": expect,
            }
        )
        if not ok:
            missing.append(path_label + f" == {expect!r}")
    return checks, missing


def run_consistency_checks(
    config: dict[str, Any], sources: dict[str, Any], source_paths: dict[str, Path]
) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    for item in config["consistency_checks"]:
        left_path = item["left"]["path"]
        right_path = item["right"]["path"]
        left_source = item["left"]["source"]
        right_source = item["right"]["source"]
        left_value = float(get_by_path(sources[left_source], left_path))
        right_value = float(get_by_path(sources[right_source], right_path))
        delta = abs(left_value - right_value)
        tolerance = float(item["tolerance"])
        results.append(
            {
                "id": item["id"],
                "left": {
                    "path": f"{relative_to_repo(source_paths[left_source])}::" + ".".join(left_path),
                    "value": left_value,
                },
                "right": {
                    "path": f"{relative_to_repo(source_paths[right_source])}::" + ".".join(right_path),
                    "value": right_value,
                },
                "abs_delta": delta,
                "tolerance": tolerance,
                "consistent": delta <= tolerance,
            }
        )
    return results


def build_context_digest(sources: dict[str, Any]) -> dict[str, Any]:
    diagrams = str(sources["diagrams"])
    return {
        "q1949": {
            "verdict": sources["q1949"].get("verdict"),
            "channel_complementarity": sources["q1949"].get("metrics", {}).get("dual", {}).get("channel_complementarity"),
        },
        "q1950": {
            "verdict": sources["q1950"].get("verdict"),
            "observer_tau": sources["q1950"].get("derived_internal_observer_params", {}).get("observer_tau"),
            "observer_feedback_gain": sources["q1950"].get("derived_internal_observer_params", {}).get("observer_feedback_gain"),
        },
        "q1951": {
            "verdict": sources["q1951"].get("verdict"),
            "mass_gain": sources["q1951"].get("mass_informational_weights", {}).get("mass_gain"),
        },
        "q1952": {
            "verdict": sources["q1952"].get("verdict"),
            "orientation_psi0": sources["q1952"].get("derived_params", {}).get("orientation_psi0"),
            "retard_phase": sources["q1952"].get("derived_params", {}).get("retard_phase"),
            "anisotropy_strength": sources["q1952"].get("derived_params", {}).get("anisotropy_strength"),
        },
        "q1953": {
            "verdict": sources["q1953"].get("verdict"),
            "tau_h": sources["q1953"].get("derived_params", {}).get("tau_h"),
            "tau_l": sources["q1953"].get("derived_params", {}).get("tau_l"),
        },
        "q1954": {
            "verdict": sources["q1954"].get("verdict"),
            "readiness": sources["q1954"].get("readiness"),
        },
        "q1955": {
            "verdict": sources["q1955"].get("verdict"),
            "repair_pass": sources["q1955"].get("repair_pass"),
        },
        "q1956": {
            "verdict": sources["q1956"].get("verdict"),
            "pass": sources["q1956"].get("pass"),
        },
        "diagrams": {
            "mentions_kernel_form": "K(d)" in diagrams and "cos(" in diagrams,
            "mentions_viscosity": "viscosity" in diagrams.lower(),
        },
    }


def build_test_results(
    config: dict[str, Any], input_values: dict[str, dict[str, Any]]
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    tolerance = float(config["trivial_selector_tolerance"])
    omega = input_values["omega"]["value"]
    retard_phase = input_values["retard_phase"]["value"]
    anis = input_values["anisotropy_strength"]["value"]
    psi0 = input_values["psi0"]["value"]
    tau_h = input_values["tau_h"]["value"]
    tau_l = input_values["tau_l"]["value"]
    speed_c = float(config["operator_rules"]["speed_of_signal_c"])

    l0 = speed_c * retard_phase / omega

    phase_iso = omega * l0 / speed_c
    lambda_iso = math.cos(phase_iso)
    matrix_a = [[lambda_iso, 0.0], [0.0, lambda_iso]]
    test_a = {
        "id": "A",
        "name": "isotropic_c_based_retardation_baseline",
        "lane": config["lane"],
        "inputs": {
            "omega": omega,
            "c": speed_c,
            "L": l0,
            "phase": phase_iso,
        },
        "expected_outcome": "a_1 = d_1 and b_1 = 0",
        "operator_form": "cos(omega * L / c) * I_2",
        "result": summarize_operator(matrix_a, tolerance),
    }

    matrix_b = [[1.0, 0.0], [0.0, 1.0]]
    test_b = {
        "id": "B",
        "name": "psi0_only_coordinate_embedding_baseline",
        "lane": config["lane"],
        "inputs": {
            "psi0": psi0,
        },
        "expected_outcome": "chart embedding only and no selector-breaking",
        "operator_form": "R(psi0) * I_2 * R(-psi0)",
        "result": summarize_operator(matrix_b, tolerance),
    }

    l_parallel_equal = l0
    l_perp_equal = l0
    matrix_c_equal = matrix_from_rotated_diagonal(
        math.cos(omega * l_parallel_equal / speed_c),
        math.cos(omega * l_perp_equal / speed_c),
        psi0,
    )
    l_parallel = l0 * (1.0 + anis)
    l_perp = l0 * (1.0 - anis)
    lambda_parallel = math.cos(omega * l_parallel / speed_c)
    lambda_perp = math.cos(omega * l_perp / speed_c)
    matrix_c_configured = matrix_from_rotated_diagonal(lambda_parallel, lambda_perp, psi0)
    test_c = {
        "id": "C",
        "name": "psi0_plus_c_based_anisotropic_retardation",
        "lane": config["lane"],
        "operator_form": "R(psi0) * diag(cos(omega * L_parallel / c), cos(omega * L_perp / c)) * R(-psi0)",
        "equal_path_control": {
            "inputs": {
                "psi0": psi0,
                "omega": omega,
                "c": speed_c,
                "L_parallel": l_parallel_equal,
                "L_perp": l_perp_equal,
            },
            "result": summarize_operator(matrix_c_equal, tolerance),
        },
        "configured_path_split": {
            "inputs": {
                "psi0": psi0,
                "omega": omega,
                "c": speed_c,
                "L0": l0,
                "L_parallel": l_parallel,
                "L_perp": l_perp,
                "anisotropy_strength": anis,
                "lambda_parallel": lambda_parallel,
                "lambda_perp": lambda_perp,
            },
            "result": summarize_operator(matrix_c_configured, tolerance),
        },
    }

    nu_parallel = 1.0 / tau_h
    nu_perp = 1.0 / tau_l
    matrix_d = matrix_from_rotated_diagonal(nu_parallel, nu_perp, psi0)
    test_d = {
        "id": "D",
        "name": "psi0_plus_viscosity_secondary_lane",
        "lane": config["lane"],
        "operator_form": "R(psi0) * diag(nu_parallel, nu_perp) * R(-psi0)",
        "inputs": {
            "psi0": psi0,
            "tau_h": tau_h,
            "tau_l": tau_l,
            "nu_parallel": nu_parallel,
            "nu_perp": nu_perp,
        },
        "result": summarize_operator(matrix_d, tolerance),
        "difference_vs_test_C_configured": summarize_operator(
            matrix_difference(matrix_d, matrix_c_configured), tolerance
        ),
    }

    primary_operator = {
        "source_test": config["selected_primary_test"],
        "lane": config["lane"],
        "operator_object": "A_1_ext(pair1)",
        **test_c["configured_path_split"]["result"],
    }
    return {
        "A": test_a,
        "B": test_b,
        "C": test_c,
        "D": test_d,
    }, primary_operator, {
        "L0": l0,
        "speed_of_signal_c": speed_c,
        "baseline_path_rule": config["operator_rules"]["baseline_path_length_rule"],
        "anisotropic_path_rule": config["operator_rules"]["anisotropic_path_length_rule"],
        "viscosity_rule": config["operator_rules"]["viscosity_rule"],
    }


def build_report(config_path: Path, report_path: Path) -> tuple[dict[str, Any], int]:
    config = load_json(config_path)
    sources, source_paths = load_sources(config)
    input_values, missing_inputs = collect_required_inputs(config, sources, source_paths)
    frontier_checks, frontier_missing = run_frontier_checks(config, sources, source_paths)
    consistency = run_consistency_checks(config, sources, source_paths)

    report: dict[str, Any] = {
        "probe_id": config["probe_id"],
        "as_of": config["as_of"],
        "lane": config["lane"],
        "pair": config["pair"],
        "basis_order": config["basis_order"],
        "config_path": relative_to_repo(config_path),
        "report_path": relative_to_repo(report_path),
        "no_false_pass": True,
        "strict_core_promotion": False,
        "context_digest": build_context_digest(sources),
        "required_inputs": input_values,
        "frozen_frontier_checks": frontier_checks,
        "source_consistency_checks": consistency,
    }

    missing_upstream_objects = sorted(set(missing_inputs + frontier_missing))
    if missing_upstream_objects:
        report.update(
            {
                "status": config["status_if_missing"],
                "selected_operator": {
                    "source_test": config["selected_primary_test"],
                    "operator_object": "A_1_ext(pair1)",
                    "classifier": "UNCOMPUTABLE",
                },
                "missing_upstream_objects": missing_upstream_objects,
                "required_next_step": "IMPLEMENT_ONE_OF_THE_LISTED_UPSTREAM_OBJECTS_BEFORE_RERUN",
                "project_conclusion": {
                    "result_type": "HARD_NEGATIVE",
                    "summary": "the repository still lacks a computable extension-lane selector operator source under the configured probe requirements",
                    "strict_core_frontier_unchanged": True,
                },
            }
        )
        return report, 1

    tests, primary_operator, operator_rules = build_test_results(config, input_values)
    primary_classifier = primary_operator["classifier"]
    if primary_classifier == "ANCHOR_IMPORTED_SPLIT":
        project_conclusion = {
            "result_type": "EXTENSION_LANE_CONCRETE_SELECTOR_SPLIT",
            "summary": "the extension lane yields a concrete selector-sector split on pair1, but the split is anchor-imported and not strict-core generated",
            "strict_core_frontier_unchanged": True,
            "required_next_step": "PHYSICAL_INTERPRETATION_OF_ANCHOR_IMPORTED_SPLIT",
        }
        required_next_step = "PHYSICAL_INTERPRETATION_OF_ANCHOR_IMPORTED_SPLIT"
    else:
        project_conclusion = {
            "result_type": "TRIVIAL_BASELINE_ONLY",
            "summary": "the configured extension probe remains selector-trivial on pair1",
            "strict_core_frontier_unchanged": True,
            "required_next_step": "SUPPLY_NONTRIVIAL_ANISOTROPIC_EXTENSION_INPUTS",
        }
        required_next_step = "SUPPLY_NONTRIVIAL_ANISOTROPIC_EXTENSION_INPUTS"

    report.update(
        {
            "status": "COMPUTED",
            "operator_rules": operator_rules,
            "tests": tests,
            "selected_operator": primary_operator,
            "missing_upstream_objects": [],
            "required_next_step": required_next_step,
            "project_conclusion": project_conclusion,
        }
    )
    return report, 0


def main() -> int:
    config_path = DEFAULT_CONFIG
    report_path = DEFAULT_REPORT
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report, exit_code = build_report(config_path, report_path)
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    summary = {
        "probe_id": report["probe_id"],
        "status": report["status"],
        "selected_operator": report.get("selected_operator", {}),
        "required_next_step": report["required_next_step"],
    }
    print(json.dumps(summary, ensure_ascii=True))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
