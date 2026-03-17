#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

AS_OF = "2026-03-17"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_PROJECTIVE_CLOSURE_F672 = GENERATED / "selector_closure_global_c_v1_projective_strict_v1.json"
IN_PSI_HESSIAN_F459 = GENERATED / "psi_hessian_diagonal_local_strict_derived_value_instantiated_v1.json"

OUT_JSON = GENERATED / "p694_current_strict_physical_computability_mass_spectrum_proxy_from_projective_selector_closure_probe.json"
OUT_SUMMARY = (
    GENERATED / "p694_current_strict_physical_computability_mass_spectrum_proxy_from_projective_selector_closure_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_finite_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def is_numeric_list_len(obj: Any, n: int) -> bool:
    return isinstance(obj, list) and len(obj) == n and all(is_finite_number(v) for v in obj)


def is_numeric_matrix(obj: Any, n: int) -> bool:
    return (
        isinstance(obj, list)
        and len(obj) == n
        and all(isinstance(row, list) and len(row) == n and all(is_finite_number(v) for v in row) for row in obj)
    )


def rayleigh_quotient(H: np.ndarray, u: np.ndarray) -> float:
    n2 = float(u @ u)
    if not (n2 > 0.0 and math.isfinite(n2)):
        raise ValueError("invalid vector norm for rayleigh quotient")
    return float((u @ (H @ u)) / n2)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_PROJECTIVE_CLOSURE_F672, IN_PSI_HESSIAN_F459]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P694",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    closure = load_json(IN_PROJECTIVE_CLOSURE_F672)
    charts = closure.get("charts")
    if not isinstance(charts, dict) or not all(isinstance(k, str) for k in charts.keys()):
        artifact = {
            "stage": "P694",
            "status": "INVALID_F672_CHARTS_SHAPE",
            "as_of": AS_OF,
            "error": "F672 must export charts as dict[str, ...]",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    expected_pairs = [f"pair{m}" for m in range(1, 6)]
    missing_pairs = [p for p in expected_pairs if p not in charts]
    if missing_pairs:
        artifact = {
            "stage": "P694",
            "status": "NOT_COMPUTABLE_MISSING_PAIR_CHARTS_IN_F672",
            "as_of": AS_OF,
            "missing_pairs": missing_pairs,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    a_refs: dict[str, Path] = {}
    missing_a_refs: list[str] = []
    for pair in expected_pairs:
        chart = charts.get(pair)
        if not isinstance(chart, dict):
            artifact = {
                "stage": "P694",
                "status": "INVALID_F672_CHART_ENTRY",
                "as_of": AS_OF,
                "error": f"F672 charts[{pair}] must be a dict",
                "no_false_pass": True,
            }
            OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            print(OUT_SUMMARY)
            return
        ref = chart.get("A_m_ref")
        if not isinstance(ref, str) or not ref:
            missing_a_refs.append(f"{pair}.A_m_ref")
            continue
        path = Path(ref)
        if not path.is_absolute():
            path = REPO / path
        if not path.exists():
            missing_a_refs.append(str(Path(ref)))
            continue
        a_refs[pair] = path

    if missing_a_refs:
        artifact = {
            "stage": "P694",
            "status": "NOT_COMPUTABLE_MISSING_PROJECTIVE_STATE_PROJECTORS",
            "as_of": AS_OF,
            "missing": missing_a_refs,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f459 = load_json(IN_PSI_HESSIAN_F459)
    out459 = f459.get("outputs") or {}
    h_psi = out459.get("h_psi")
    evals = out459.get("eigenvalues_ascending")
    ecols = out459.get("eigenvectors_columns")

    if not is_numeric_matrix(h_psi, 12):
        artifact = {
            "stage": "P694",
            "status": "INVALID_F459_H_PSI_SHAPE",
            "as_of": AS_OF,
            "error": "F459 outputs.h_psi must be a 12x12 numeric matrix",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return
    if not (isinstance(evals, list) and len(evals) == 12 and all(is_finite_number(x) for x in evals)):
        artifact = {
            "stage": "P694",
            "status": "INVALID_F459_EIGENVALUES_SHAPE",
            "as_of": AS_OF,
            "error": "F459 outputs.eigenvalues_ascending must be length-12 numeric list",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return
    if not (isinstance(ecols, list) and len(ecols) == 12 and all(is_numeric_list_len(col, 12) for col in ecols)):
        artifact = {
            "stage": "P694",
            "status": "INVALID_F459_EIGENVECTORS_SHAPE",
            "as_of": AS_OF,
            "error": "F459 outputs.eigenvectors_columns must be list of 12 length-12 numeric lists",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    H = np.array(h_psi, dtype=float)
    sym_residual_inf = float(np.max(np.abs(H - H.T)))
    Hs = 0.5 * (H + H.T)

    u_by_pair: dict[str, np.ndarray] = {}
    u_norm_by_pair: dict[str, float] = {}
    theta_by_pair: dict[str, float | None] = {}
    for pair, a_path in a_refs.items():
        a = load_json(a_path)
        data = a.get("data") or {}
        m = int(pair.replace("pair", ""))
        u_key = f"u_{m}"
        u = data.get(u_key)
        if not is_numeric_list_len(u, 12):
            artifact = {
                "stage": "P694",
                "status": "INVALID_A_M_U_VECTOR_SHAPE",
                "as_of": AS_OF,
                "error": f"{a_path.relative_to(REPO)} must contain data.{u_key} as length-12 numeric list",
                "no_false_pass": True,
            }
            OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            print(OUT_SUMMARY)
            return
        uv = np.array([float(x) for x in u], dtype=float)
        u_by_pair[pair] = uv
        u_norm_by_pair[pair] = float(np.linalg.norm(uv))

        t = data.get(f"theta_{m}_exported")
        theta_by_pair[pair] = float(t) if is_finite_number(t) else None

    m2_by_pair: dict[str, float] = {}
    for pair in expected_pairs:
        m2_by_pair[pair] = float(rayleigh_quotient(Hs, u_by_pair[pair]))

    # Eigenmode overlaps against the selector rays u_m (sign-gauge-invariant via abs^2).
    V = np.column_stack([np.array([float(x) for x in col], dtype=float) for col in ecols])
    v_norms = np.linalg.norm(V, axis=0)
    v_norm_max_abs_diff_from_1 = float(np.max(np.abs(v_norms - 1.0)))

    per_eigenmode: list[dict[str, Any]] = []
    for j in range(12):
        v = V[:, j]
        v_n2 = float(v @ v)
        overlap_abs2: dict[str, float] = {}
        for pair in expected_pairs:
            u = u_by_pair[pair]
            u_n2 = float(u @ u)
            dot = float(u @ v)
            overlap_abs2[pair] = float((dot * dot) / (u_n2 * v_n2))
        dom_pair = max(overlap_abs2.items(), key=lambda kv: kv[1])
        per_eigenmode.append(
            {
                "j": int(j),
                "lambda": float(evals[j]),
                "dominant_pair": {"pair": dom_pair[0], "overlap_abs2": float(dom_pair[1])},
                "overlap_abs2_by_pair": overlap_abs2,
            }
        )

    m2_sorted = sorted(((p, float(v)) for p, v in m2_by_pair.items()), key=lambda kv: kv[1])

    status = "PASS_PHYSICAL_MASS_SPECTRUM_PROXY_COMPUTABLE_FROM_PROJECTIVE_SELECTOR_CLOSURE"
    artifact = {
        "stage": "P694",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "export_a_probe_level_physical_computability_witness__compute_a_first_quadratic_mass_spectrum_proxy_from_the_exported_global_projective_selector_closure_rays_on_C_v1_and_the_exported_diagonal_local_H_psi_value_instantiation__sign_gauge_invariant__no_false_pass",
        "inputs": {
            "projective_selector_closure_ref": str(IN_PROJECTIVE_CLOSURE_F672.relative_to(REPO)),
            "psi_hessian_ref": str(IN_PSI_HESSIAN_F459.relative_to(REPO)),
            "orientation_projector_operators_by_pair": {p: str(a_refs[p].relative_to(REPO)) for p in expected_pairs},
        },
        "audits": {
            "h_psi_symmetry_residual_inf_norm": sym_residual_inf,
            "u_norm_by_pair": {p: float(u_norm_by_pair[p]) for p in expected_pairs},
            "f459_eigenvector_norm_max_abs_diff_from_1": v_norm_max_abs_diff_from_1,
        },
        "results": {
            "theta_by_pair_if_present": theta_by_pair,
            "mass_proxy_m2_by_pair_rayleigh_on_H_psi": m2_by_pair,
            "mass_proxy_m2_sorted_ascending": [{"pair": p, "m2": float(v)} for p, v in m2_sorted],
            "h_psi_eigenvalues_ascending": [float(x) for x in evals],
            "per_eigenmode_overlap_against_selector_rays": per_eigenmode,
        },
        "hard_limits": [
            "Probe-only operational computability witness; does not claim Standard Model identification of any computed mass proxy.",
            "Uses a specific diagonal/local strict-derived H_psi value instantiation (F459); does not promote that instantiation to host matching or global physics.",
            "Uses projective/ray selector data only (u and -u identified); does not claim any directed/sign-sensitive physical orientation datum in strict core.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P694",
        "status": status,
        "mass_proxy_m2_by_pair_rayleigh_on_H_psi": m2_by_pair,
        "mass_proxy_m2_sorted_ascending": [{"pair": p, "m2": float(v)} for p, v in m2_sorted],
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

