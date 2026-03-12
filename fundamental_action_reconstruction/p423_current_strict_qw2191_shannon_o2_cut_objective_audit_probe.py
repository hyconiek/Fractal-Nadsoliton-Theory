#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_QW2190 = REPO / "report_qw2190_kernel_mode_representation_emergence_gate.json"

OUT_JSON = GENERATED / "p423_current_strict_qw2191_shannon_o2_cut_objective_audit_probe.json"
OUT_SUMMARY = GENERATED / "p423_current_strict_qw2191_shannon_o2_cut_objective_audit_probe_summary.json"


def load_json_path(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def real_fourier_basis(n: int) -> dict[str, np.ndarray]:
    x = np.arange(n, dtype=float)
    basis: dict[str, np.ndarray] = {"e0": np.ones(n, dtype=float) / math.sqrt(n)}
    for m in range(1, n // 2):
        basis[f"c{m}"] = math.sqrt(2.0 / n) * np.cos(2.0 * math.pi * m * x / n)
        basis[f"s{m}"] = math.sqrt(2.0 / n) * np.sin(2.0 * math.pi * m * x / n)
    if n % 2 == 0:
        basis[f"e{n//2}"] = ((-1.0) ** x) / math.sqrt(n)
    return basis


def rotate_pair(c: np.ndarray, s: np.ndarray, theta: float) -> tuple[np.ndarray, np.ndarray]:
    ct = math.cos(theta)
    st = math.sin(theta)
    return ct * c + st * s, -st * c + ct * s


def prob_from_vector_sq(u: np.ndarray) -> np.ndarray:
    p = np.array(u, dtype=float) ** 2
    z = float(np.sum(p))
    if z <= 0:
        raise ValueError("vector square sum must be positive")
    return p / z


def shannon_H(p: np.ndarray) -> float:
    eps = 1e-18
    p2 = np.clip(np.array(p, dtype=float), eps, 1.0)
    p2 = p2 / float(np.sum(p2))
    return float(-np.sum(p2 * np.log(p2)))


def kl(p: np.ndarray, q: np.ndarray) -> float:
    eps = 1e-18
    p2 = np.clip(np.array(p, dtype=float), eps, 1.0)
    q2 = np.clip(np.array(q, dtype=float), eps, 1.0)
    p2 = p2 / float(np.sum(p2))
    q2 = q2 / float(np.sum(q2))
    return float(np.sum(p2 * np.log(p2 / q2)))


def sym_kl(p: np.ndarray, q: np.ndarray) -> float:
    return float(kl(p, q) + kl(q, p))


@dataclass(frozen=True)
class ObjectiveSpec:
    id: str
    description: str
    fn: Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray], float]


def cluster_near_minima(thetas: np.ndarray, values: np.ndarray) -> dict[str, Any]:
    v = np.array(values, dtype=float)
    th = np.array(thetas, dtype=float)
    if v.shape != th.shape:
        raise ValueError("thetas/values shape mismatch")

    vmin = float(np.min(v))
    vmax = float(np.max(v))
    vrng = float(vmax - vmin)
    tol = float(max(1e-12, 1e-6 * max(1.0, vrng)))

    near = np.where(v <= (vmin + tol))[0].tolist()
    if not near:
        raise RuntimeError("no near-min indices found (unexpected)")

    # Merge wrap-around cluster if needed (theta grid is periodic).
    if 0 in near and (len(v) - 1) in near:
        # rotate indices so cluster starts after the wrap.
        k = 0
        while k < len(v) and k in near:
            k += 1
        near = [i for i in near if i >= k] + [i for i in near if i < k]

    clusters: list[list[int]] = []
    cur: list[int] = []
    for idx in near:
        if not cur or idx == (cur[-1] + 1):
            cur.append(idx)
        else:
            clusters.append(cur)
            cur = [idx]
    if cur:
        clusters.append(cur)

    centers = [float(th[int(round(float(np.mean(c))))]) for c in clusters]
    return {
        "min_value": vmin,
        "max_value": vmax,
        "range": vrng,
        "near_min_tol": tol,
        "near_min_index_count": int(len(near)),
        "cluster_count": int(len(clusters)),
        "cluster_centers_theta": centers,
    }


def main() -> None:
    GENERATED.mkdir(exist_ok=True)
    q2190 = load_json_path(IN_QW2190)
    n = int(q2190["mode_mapping"]["n_octaves"])

    fb = real_fourier_basis(n)
    c1, s1 = fb["c1"], fb["s1"]
    c2, s2 = fb["c2"], fb["s2"]

    # Dense deterministic grid (audit only; not used as a selector itself).
    grid_n = 7201
    thetas = np.linspace(0.0, 2.0 * math.pi, grid_n)

    objectives: list[ObjectiveSpec] = [
        ObjectiveSpec(
            id="H_avg_u1_v1_u2_v2",
            description="Average Shannon entropy of squared site-amplitudes for (u1,v1,u2,v2) on the 12-octave ring site basis",
            fn=lambda u1, v1, u2, v2: float(
                0.25
                * (
                    shannon_H(prob_from_vector_sq(u1))
                    + shannon_H(prob_from_vector_sq(v1))
                    + shannon_H(prob_from_vector_sq(u2))
                    + shannon_H(prob_from_vector_sq(v2))
                )
            ),
        ),
        ObjectiveSpec(
            id="symKL_u1_u2",
            description="Symmetric KL divergence between squared site-amplitude distributions of u1 and u2",
            fn=lambda u1, _v1, u2, _v2: float(
                sym_kl(prob_from_vector_sq(u1), prob_from_vector_sq(u2))
            ),
        ),
        ObjectiveSpec(
            id="symKL_u1_v2",
            description="Symmetric KL divergence between squared site-amplitude distributions of u1 and v2",
            fn=lambda u1, _v1, _u2, v2: float(
                sym_kl(prob_from_vector_sq(u1), prob_from_vector_sq(v2))
            ),
        ),
        ObjectiveSpec(
            id="H_mix_u1_u2",
            description="Shannon entropy of the equal-weight mixture (u1^2 + u2^2)/2 over ring sites",
            fn=lambda u1, _v1, u2, _v2: float(
                shannon_H(0.5 * (prob_from_vector_sq(u1) + prob_from_vector_sq(u2)))
            ),
        ),
    ]

    rows: list[dict[str, Any]] = []
    summaries: dict[str, Any] = {}

    for obj in objectives:
        vals: list[float] = []
        for th in thetas.tolist():
            u1, v1 = rotate_pair(c1, s1, float(th))
            u2, v2 = rotate_pair(c2, s2, float(th))
            vals.append(float(obj.fn(u1, v1, u2, v2)))

        arr = np.array(vals, dtype=float)
        summary = cluster_near_minima(thetas, arr)
        summaries[obj.id] = {
            "description": obj.description,
            "summary": summary,
        }

        rows.append(
            {
                "objective_id": obj.id,
                "min_value": float(summary["min_value"]),
                "cluster_count": int(summary["cluster_count"]),
                "cluster_centers_theta": [float(x) for x in summary["cluster_centers_theta"]],
            }
        )

    # A strict-core O(2)-cut ingredient would require cluster_count == 1 (mod declared residual symmetry).
    any_unique = bool(any(int(r["cluster_count"]) == 1 for r in rows))

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P423",
        "goal": "audit_naive_shannon_objective_candidates_on_qw2191_o2_family_for_unique_o2_cut",
        "inputs": {
            "qw2190_kernel_mode_report": str(IN_QW2190.relative_to(REPO)),
        },
        "grid": {"theta_domain": "[0,2pi]", "grid_n": int(grid_n)},
        "objectives": {obj.id: obj.description for obj in objectives},
        "objective_minima_summary": summaries,
        "objective_rows_compact": rows,
        "verdict": {
            "any_objective_unique_on_grid": bool(any_unique),
            "statement": (
                "On the audited dense theta grid, the tested naive Shannon-type objective candidates "
                "(site-amplitude entropies and simple KL-based divergences) do not yield a unique global minimizer "
                "on the QW-2191 O(2) family. Therefore these naive Shannon objective forms do not, by themselves, "
                "supply a strict-core canonical O(2)-cut ingredient."
                if not any_unique
                else (
                    "At least one tested objective has a single near-minimum cluster on the audited theta grid. "
                    "This probe does not promote that into strict-core evidence; it only records the grid-audit "
                    "outcome and leaves open whether a fully typed, quotient-safe uniqueness theorem exists."
                )
            ),
        },
        "hard_limits": [
            "Audit only: does not export a strict-core selector ingredient.",
            "Does not discharge QW-2191.",
            "Does not claim a strict Shannon symmetry-breaking law exists.",
            "Does not discharge T165 (Shannon symmetry-breaking selector ingredient target).",
            "Does not claim strict-core theta export or ToE closure.",
        ],
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    summary_out = {
        "stage": report["stage"],
        "any_objective_unique_on_grid": bool(any_unique),
        "objective_rows_compact": rows,
        "no_false_pass": True,
    }
    OUT_SUMMARY.write_text(json.dumps(summary_out, ensure_ascii=False, indent=2), encoding="utf-8")

    print(json.dumps({"stage": "P423", "any_unique": bool(any_unique)}, ensure_ascii=False))


if __name__ == "__main__":
    main()

