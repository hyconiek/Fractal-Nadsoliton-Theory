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
IN_F92 = GENERATED / "f92_first_preobserver_to_emergent_observer_coarse_graining_packet_summary.json"
IN_F93 = GENERATED / "f93_first_emergent_observer_limit_readout_operator_packet_summary.json"
IN_F94 = GENERATED / "f94_first_actual_emergent_observer_construction_candidate_packet_summary.json"
IN_F95 = GENERATED / "f95_first_actual_emergent_observer_realization_map_packet_summary.json"
IN_F96 = GENERATED / "f96_first_emergent_observer_self_consistency_operator_packet_summary.json"
IN_F97 = GENERATED / "f97_first_emergent_observer_fixed_point_reduction_packet_summary.json"
IN_F98 = GENERATED / "f98_first_emergent_observer_fixed_point_object_candidate_packet_summary.json"
IN_F99 = GENERATED / "f99_first_emergent_observer_closure_candidate_packet_summary.json"

OUT_JSON = (
    GENERATED
    / "p699_current_strict_projective_emergent_observer_chain_from_global_projective_selector_closure_output_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p699_current_strict_projective_emergent_observer_chain_from_global_projective_selector_closure_output_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_finite_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def is_numeric_matrix(obj: Any, rows: int, cols: int) -> bool:
    return (
        isinstance(obj, list)
        and len(obj) == rows
        and all(isinstance(row, list) and len(row) == cols and all(is_finite_number(v) for v in row) for row in obj)
    )


def induced_projector_pushforward(matrix: np.ndarray, proj: np.ndarray) -> np.ndarray:
    return matrix @ proj @ matrix.T


def projector_audit(proj: np.ndarray) -> dict[str, float | list[float]]:
    sym_residual_inf = float(np.max(np.abs(proj - proj.T)))
    evals = np.linalg.eigvalsh(0.5 * (proj + proj.T))
    trace = float(np.trace(proj))
    fro = float(np.linalg.norm(proj))
    return {
        "trace": trace,
        "fro_norm": fro,
        "symmetry_residual_inf_norm": sym_residual_inf,
        "eigenvalues_symmetrized": [float(evals[0]), float(evals[-1])] if evals.size else [],
    }


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        IN_PROJECTIVE_CLOSURE_F672,
        IN_F92,
        IN_F93,
        IN_F94,
        IN_F95,
        IN_F96,
        IN_F97,
        IN_F98,
        IN_F99,
    ]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P699",
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
    out_obs = closure.get("output_observable")
    if not isinstance(out_obs, dict) or out_obs.get("basis") != ["o_+", "o_-"]:
        artifact = {
            "stage": "P699",
            "status": "INVALID_F672_OUTPUT_OBSERVABLE",
            "as_of": AS_OF,
            "error": "F672 must export output_observable.basis = ['o_+','o_-']",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    P_out = out_obs.get("output_projector_matrix_in_o_plus_o_minus")
    if not is_numeric_matrix(P_out, 2, 2):
        artifact = {
            "stage": "P699",
            "status": "INVALID_F672_OUTPUT_PROJECTOR_SHAPE",
            "as_of": AS_OF,
            "error": "F672 output_observable.output_projector_matrix_in_o_plus_o_minus must be a 2x2 numeric matrix",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f92 = load_json(IN_F92)
    f93 = load_json(IN_F93)
    f94 = load_json(IN_F94)
    f95 = load_json(IN_F95)
    f96 = load_json(IN_F96)
    f97 = load_json(IN_F97)
    f98 = load_json(IN_F98)
    f99 = load_json(IN_F99)

    C = ((f92.get("coarse_graining_operator") or {}).get("matrix")) if isinstance(f92, dict) else None
    L = ((f93.get("observer_limit_readout_operator") or {}).get("matrix")) if isinstance(f93, dict) else None
    G = ((f94.get("observer_construction_candidate_operator") or {}).get("matrix")) if isinstance(f94, dict) else None
    H = ((f95.get("observer_realization_map") or {}).get("matrix")) if isinstance(f95, dict) else None
    J = ((f96.get("observer_self_consistency_operator") or {}).get("matrix")) if isinstance(f96, dict) else None
    K = ((f97.get("observer_fixed_point_reduction_operator") or {}).get("matrix")) if isinstance(f97, dict) else None
    M = ((f98.get("observer_fixed_point_object_map") or {}).get("matrix")) if isinstance(f98, dict) else None
    N = ((f99.get("observer_closure_candidate_map") or {}).get("matrix")) if isinstance(f99, dict) else None

    if not is_numeric_matrix(C, 2, 2):
        raise SystemExit("F92 coarse_graining_operator.matrix must be 2x2")
    if not is_numeric_matrix(L, 2, 2):
        raise SystemExit("F93 observer_limit_readout_operator.matrix must be 2x2")
    if not is_numeric_matrix(G, 2, 2):
        raise SystemExit("F94 observer_construction_candidate_operator.matrix must be 2x2")
    if not is_numeric_matrix(H, 2, 2):
        raise SystemExit("F95 observer_realization_map.matrix must be 2x2")
    if not is_numeric_matrix(J, 2, 2):
        raise SystemExit("F96 observer_self_consistency_operator.matrix must be 2x2")
    if not is_numeric_matrix(K, 1, 2):
        raise SystemExit("F97 observer_fixed_point_reduction_operator.matrix must be 1x2")
    if not is_numeric_matrix(M, 1, 1):
        raise SystemExit("F98 observer_fixed_point_object_map.matrix must be 1x1")
    if not is_numeric_matrix(N, 1, 1):
        raise SystemExit("F99 observer_closure_candidate_map.matrix must be 1x1")

    Cn = np.array(C, dtype=float)
    Ln = np.array(L, dtype=float)
    Gn = np.array(G, dtype=float)
    Hn = np.array(H, dtype=float)
    Jn = np.array(J, dtype=float)
    Kn = np.array(K, dtype=float)
    Mn = np.array(M, dtype=float)
    Nn = np.array(N, dtype=float)

    P0 = np.array(P_out, dtype=float)
    Py = induced_projector_pushforward(Cn, P0)
    Pz = induced_projector_pushforward(Ln, Py)
    Pw = induced_projector_pushforward(Gn, Pz)
    Px = induced_projector_pushforward(Hn, Pw)
    Pu = induced_projector_pushforward(Jn, Px)
    Pf = induced_projector_pushforward(Kn, Pu)
    Pp = induced_projector_pushforward(Mn, Pf)
    Pc = induced_projector_pushforward(Nn, Pp)

    tol = 1e-12

    def _diag2(a: np.ndarray) -> tuple[float, float]:
        return float(a[0, 0]), float(a[1, 1])

    z_commit, z_residual = _diag2(Pz)
    w_commit, w_residual = _diag2(Pw)
    x_commit, x_residual = _diag2(Px)
    u_commit, u_residual = _diag2(Pu)
    c_closure = float(Pc[0, 0])

    ok = True
    ok = ok and (z_commit > 0.0) and (abs(z_residual) <= tol)
    ok = ok and (w_commit > 0.0) and (abs(w_residual) <= tol)
    ok = ok and (x_commit > 0.0) and (abs(x_residual) <= tol)
    ok = ok and (u_commit > 0.0) and (abs(u_residual) <= tol)
    ok = ok and (c_closure > 0.0)

    status = "PASS_PROJECTIVE_EMERGENT_OBSERVER_CHAIN_COMPUTABLE_FROM_PROJECTIVE_SELECTOR_CLOSURE_OUTPUT"
    if not ok:
        status = "FAIL_PROJECTIVE_EMERGENT_OBSERVER_CHAIN_NOT_COMPUTABLE_OR_RESIDUAL_NONVANISHING"

    artifact = {
        "stage": "P699",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "compute_projective_emergent_observer_chain_from_global_projective_selector_closure_output_projector__treat_directed_sign_as_gauge__no_false_pass",
        "inputs": {
            "projective_selector_closure_ref": str(IN_PROJECTIVE_CLOSURE_F672.relative_to(REPO)),
            "F92": str(IN_F92.relative_to(REPO)),
            "F93": str(IN_F93.relative_to(REPO)),
            "F94": str(IN_F94.relative_to(REPO)),
            "F95": str(IN_F95.relative_to(REPO)),
            "F96": str(IN_F96.relative_to(REPO)),
            "F97": str(IN_F97.relative_to(REPO)),
            "F98": str(IN_F98.relative_to(REPO)),
            "F99": str(IN_F99.relative_to(REPO)),
        },
        "operators": {
            "C_obs_limit_preLM_v1": {"matrix": C},
            "L_obs_limit_preLM_v1": {"matrix": L},
            "G_obs_candidate_preLM_v1": {"matrix": G},
            "H_obs_realization_preLM_v1": {"matrix": H},
            "J_obs_self_consistency_preLM_v1": {"matrix": J},
            "K_obs_fixed_point_preLM_v1": {"matrix": K},
            "M_obs_fixed_object_preLM_v1": {"matrix": M},
            "N_obs_closure_candidate_preLM_v1": {"matrix": N},
        },
        "projectors": {
            "Q_out(o_+,o_-)": {"matrix": P_out, "audit": projector_audit(P0)},
            "Y_obs_limit(y_bias,y_total)": {"matrix": Py.tolist(), "audit": projector_audit(Py)},
            "Z_obs_limit(z_commit,z_residual)": {"matrix": Pz.tolist(), "audit": projector_audit(Pz)},
            "W_obs_candidate(w_commit,w_residual)": {"matrix": Pw.tolist(), "audit": projector_audit(Pw)},
            "X_obs_real(x_commit,x_residual)": {"matrix": Px.tolist(), "audit": projector_audit(Px)},
            "U_obs_cons(u_commit,u_residual)": {"matrix": Pu.tolist(), "audit": projector_audit(Pu)},
            "F_obs_fix(f_commit)": {"matrix": Pf.tolist(), "audit": projector_audit(Pf)},
            "P_obs_fix_obj(p_fix)": {"matrix": Pp.tolist(), "audit": projector_audit(Pp)},
            "C_obs_closure(c_closure)": {"matrix": Pc.tolist(), "audit": projector_audit(Pc)},
        },
        "results": {
            "z_commit": z_commit,
            "z_residual": z_residual,
            "w_commit": w_commit,
            "w_residual": w_residual,
            "x_commit": x_commit,
            "x_residual": x_residual,
            "u_commit": u_commit,
            "u_residual": u_residual,
            "c_closure": c_closure,
            "tolerance_abs": tol,
        },
        "hard_limits": [
            "Projective/gauge-safe downstream chain only; does not claim any directed/sign-sensitive physical orientation datum in strict core.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim strict-core or global ToE closure.",
            "Does not claim actual emergent observer closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P699",
        "status": status,
        "z_commit": z_commit,
        "z_residual": z_residual,
        "u_commit": u_commit,
        "u_residual": u_residual,
        "c_closure": c_closure,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

