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
IN_R14 = GENERATED / "r14_explicit_frozen_kernel_specialization_packet_for_host_matching_route.json"

T165_THETA_FIX = GENERATED / "theta_fix_pair1_o2_cut_ord_reference_v1.json"

OUT_JSON = GENERATED / "p439_current_strict_qw2191_weighted_kl_reference_objective_o2_cut_audit_probe.json"
OUT_SUMMARY = GENERATED / "p439_current_strict_qw2191_weighted_kl_reference_objective_o2_cut_audit_probe_summary.json"


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


def normalize_prob(v: np.ndarray) -> np.ndarray:
    vv = np.array(v, dtype=float)
    if np.any(vv < 0):
        raise ValueError("reference weights must be nonnegative")
    z = float(np.sum(vv))
    if z <= 0:
        raise ValueError("reference weights must sum to positive")
    return vv / z


def kl(p: np.ndarray, q: np.ndarray) -> float:
    eps = 1e-18
    p2 = np.clip(np.array(p, dtype=float), eps, 1.0)
    q2 = np.clip(np.array(q, dtype=float), eps, 1.0)
    p2 = p2 / float(np.sum(p2))
    q2 = q2 / float(np.sum(q2))
    return float(np.sum(p2 * np.log(p2 / q2)))


@dataclass(frozen=True)
class RefSpec:
    id: str
    description: str
    weights: np.ndarray


@dataclass(frozen=True)
class ObjectiveSpec:
    id: str
    description: str
    fn: Callable[[np.ndarray, np.ndarray, np.ndarray], float]


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


def extract_ktotal_from_r14(r14: dict[str, Any]) -> np.ndarray:
    rows = r14.get("specialization_rows")
    if not isinstance(rows, list) or not rows:
        raise ValueError("R14 specialization_rows missing or empty")
    n = len(rows)
    m = np.zeros((n, n), dtype=float)
    for i, row in enumerate(rows):
        if not isinstance(row, list) or len(row) != n:
            raise ValueError("R14 specialization_rows not square")
        for j, entry in enumerate(row):
            if not isinstance(entry, dict) or "specialized_value" not in entry:
                raise ValueError("R14 specialization entry missing specialized_value")
            m[i, j] = float(entry["specialized_value"])
    return m


def centers_form_z2_pair(centers: list[float], tol: float = 1e-3) -> bool:
    if len(centers) != 2:
        return False
    a, b = float(min(centers)), float(max(centers))
    diff = (b - a) % (2.0 * math.pi)
    if diff > math.pi:
        diff = 2.0 * math.pi - diff
    return abs(diff - math.pi) <= tol


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    q2190 = load_json_path(IN_QW2190)
    n = int(q2190["mode_mapping"]["n_octaves"])
    if n != 12:
        raise ValueError(f"expected n=12 from QW-2190 scaffold, got n={n}")

    r14 = load_json_path(IN_R14)
    ktotal = extract_ktotal_from_r14(r14)
    if ktotal.shape != (n, n):
        raise ValueError("unexpected K_total shape from R14")

    fb = real_fourier_basis(n)
    c1, s1 = fb["c1"], fb["s1"]

    alpha_geo = 4.0 * math.log(2.0)
    sigma_int = -1.0

    # Build a small family of explicit non-uniform reference distributions r(x).
    # NOTE: These references are *not* claimed to be strict-core admissible selector ingredients.
    # They are audited only as objective-shape candidates for the T165 lane (N464 escape hatch).
    x = np.arange(n, dtype=float)
    d_cyc = np.minimum(x, n - x)  # cyclic distance from site 0
    elem_order = np.array([1 if i == 0 else n // math.gcd(int(i), n) for i in range(n)], dtype=float)

    refs: list[RefSpec] = [
        RefSpec(
            id="r_uniform",
            description="Uniform reference distribution on ring sites (control; expected to inherit N463/N464 periodicity)",
            weights=np.ones(n, dtype=float),
        ),
        RefSpec(
            id="r_1_plus_abs_Ktotal_row0",
            description="Reference weights w_x := 1 + |K_total[0,x]| from R14 (marked-site non-uniform datum; normalized to prob)",
            weights=1.0 + np.abs(ktotal[0, :]),
        ),
        RefSpec(
            id="r_exp_minus_alpha_geo_dcyc",
            description="Reference weights w_x := exp(-alpha_geo * d_cyc(0,x)) with alpha_geo=4 ln 2 (marked-site distance decay; normalized)",
            weights=np.exp(-alpha_geo * d_cyc),
        ),
        RefSpec(
            id="r_exp_minus_alpha_geo_directed_d",
            description="Reference weights w_x := exp(-alpha_geo * d_dir(0,x)) where d_dir(0,x)=x on Z_12 (marked-site + marked-direction; normalized)",
            weights=np.exp(-alpha_geo * x),
        ),
        RefSpec(
            id="r_exp_sigma_int_alpha_geo_Ktotal_row0",
            description="Reference weights w_x := exp(sigma_int * alpha_geo * K_total[0,x]) with sigma_int=-1, alpha_geo=4 ln 2 (kernel-shaped; normalized)",
            weights=np.exp(sigma_int * alpha_geo * ktotal[0, :]),
        ),
        RefSpec(
            id="r_exp_minus_alpha_geo_element_order",
            description=(
                "Reference weights w_x := exp(-alpha_geo * ord_Z12(x)) where ord_Z12(x)=order of element x in Z_12 "
                "(group-structure-derived non-uniform datum; no marked direction; normalized)"
            ),
            weights=np.exp(-alpha_geo * elem_order),
        ),
        RefSpec(
            id="r_element_order_power_minus_alpha_geo",
            description=(
                "Reference weights w_x := ord_Z12(x)^(-alpha_geo) with alpha_geo=4 ln 2 "
                "(group-structure-derived power-law; no marked direction; normalized)"
            ),
            weights=np.exp(-alpha_geo * np.log(elem_order)),
        ),
    ]

    refs_norm = {r.id: normalize_prob(r.weights) for r in refs}

    objectives: list[ObjectiveSpec] = [
        ObjectiveSpec(
            id="KL_u1_to_r",
            description="J(theta)=KL(p_theta || r) where p_theta=u1(theta)^2 is the squared site-amplitude distribution",
            fn=lambda p, _q, r: float(kl(p, r)),
        ),
        ObjectiveSpec(
            id="KL_avg_u1_v1_to_r",
            description="J(theta)=0.5*( KL(p_theta||r) + KL(q_theta||r) ) where q_theta=v1(theta)^2 (orthogonal partner)",
            fn=lambda p, q, r: float(0.5 * (kl(p, r) + kl(q, r))),
        ),
    ]

    grid_n = 7201
    thetas = np.linspace(0.0, 2.0 * math.pi, grid_n)

    objective_table: dict[str, Any] = {}
    any_unique_cluster = False
    any_z2_unique = False

    for ref in refs:
        rprob = refs_norm[ref.id]
        ref_rows: list[dict[str, Any]] = []
        ref_details: dict[str, Any] = {}
        for obj in objectives:
            vals: list[float] = []
            for th in thetas.tolist():
                u1, v1 = rotate_pair(c1, s1, float(th))
                p = prob_from_vector_sq(u1)
                q = prob_from_vector_sq(v1)
                vals.append(float(obj.fn(p, q, rprob)))

            arr = np.array(vals, dtype=float)
            summary = cluster_near_minima(thetas, arr)
            cc = int(summary["cluster_count"])
            centers = [float(x) for x in summary["cluster_centers_theta"]]
            z2_unique = bool(cc == 2 and centers_form_z2_pair(centers))

            any_unique_cluster = any_unique_cluster or (cc == 1)
            any_z2_unique = any_z2_unique or z2_unique
            ref_details[obj.id] = {
                "description": obj.description,
                "summary": summary,
                "z2_unique_on_grid": z2_unique,
            }
            ref_rows.append(
                {
                    "objective_id": obj.id,
                    "min_value": float(summary["min_value"]),
                    "cluster_count": cc,
                    "cluster_centers_theta": centers,
                    "z2_unique_on_grid": z2_unique,
                }
            )

        objective_table[ref.id] = {
            "reference_description": ref.description,
            "reference_prob": [float(v) for v in rprob.tolist()],
            "objective_minima_summary": ref_details,
            "objective_rows_compact": ref_rows,
        }

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P439",
        "goal": "audit_weighted_KL_to_nonuniform_reference_objectives_on_qw2191_o2_family_for_potential_unique_o2_cut",
        "inputs": {
            "qw2190_kernel_mode_report": str(IN_QW2190.relative_to(REPO)),
            "r14_ktotal_specialization_packet": str(IN_R14.relative_to(REPO)),
            "alpha_geo_used": alpha_geo,
            "sigma_int_used": sigma_int,
        },
        "grid": {"theta_domain": "[0,2pi]", "grid_n": int(grid_n)},
        "references_tested": {r.id: r.description for r in refs},
        "objectives_tested": {o.id: o.description for o in objectives},
        "results_by_reference": objective_table,
        "verdict": {
            "any_objective_unique_cluster_on_grid": bool(any_unique_cluster),
            "any_objective_z2_unique_on_grid": bool(any_z2_unique),
            "statement": (
                "At least one weighted KL-to-reference objective yields exactly two near-minimum clusters separated by π on the audited grid "
                "(i.e. Z2-unique up to the unavoidable sign flip). This remains a probe-level observation only; it is not promoted into a strict-core "
                "O(2)-cut ingredient without (1) an explicit admissibility story for the non-translation-invariant reference datum and "
                "(2) a theorem-level uniqueness proof for the exact objective/reference pair under discussion. "
                "Note: T165 is discharged by a separate strict theorem/packet (F446/N480) for the element-order reference cross-entropy objective; "
                "this probe continues to audit alternative KL-shaped objectives/references without promotion."
                if any_z2_unique
                else (
                    "On the audited theta grid, none of the tested weighted KL-to-reference objectives yield a Z2-unique cut (two minima separated by π) "
                    "and none yield a single unique minimizer cluster. Therefore these reference-shape candidates (as tested) still do not supply a viable "
                    "strict-core O(2)-cut route for T165."
                )
            ),
        },
        "hard_limits": [
            "Audit only: does not export a strict-core selector ingredient.",
            "Does not discharge QW-2191.",
            "Does not decide T166/T168/T169.",
            "Does not claim any marked-site / reference datum is already strict-admissible as physics.",
            "Does not claim strict-core theta export or ToE closure.",
        ],
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    summary_out = {
        "stage": report["stage"],
        "any_objective_unique_cluster_on_grid": bool(any_unique_cluster),
        "any_objective_z2_unique_on_grid": bool(any_z2_unique),
        "t165_selector_ingredient_exported": bool(T165_THETA_FIX.exists()),
        "references_tested": {r.id: r.description for r in refs},
        "objectives_tested": {o.id: o.description for o in objectives},
        "no_false_pass": True,
    }
    OUT_SUMMARY.write_text(json.dumps(summary_out, ensure_ascii=False, indent=2), encoding="utf-8")

    print(
        json.dumps(
            {
                "stage": "P439",
                "any_unique_cluster": bool(any_unique_cluster),
                "any_z2_unique": bool(any_z2_unique),
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
