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

OUT_JSON = (
    GENERATED
    / "p696_current_strict_physical_computability_selector_aligned_channel_spectrum_proxy_from_projective_selector_closure_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p696_current_strict_physical_computability_selector_aligned_channel_spectrum_proxy_from_projective_selector_closure_probe_summary.json"
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


def fourier_pair_basis(n: int, m: int) -> tuple[np.ndarray, np.ndarray]:
    # Orthonormal real Fourier pair basis on Z_n: c_m, s_m with scale sqrt(2/n).
    scale = math.sqrt(2.0 / float(n))
    c = np.array([scale * math.cos(2.0 * math.pi * float(m) * float(k) / float(n)) for k in range(n)], dtype=float)
    s = np.array([scale * math.sin(2.0 * math.pi * float(m) * float(k) / float(n)) for k in range(n)], dtype=float)
    return c, s


def fourier_e0(n: int) -> np.ndarray:
    scale = 1.0 / math.sqrt(float(n))
    return np.array([scale for _ in range(n)], dtype=float)


def fourier_e6_on_n12() -> np.ndarray:
    n = 12
    scale = 1.0 / math.sqrt(float(n))
    return np.array([scale * (1.0 if (k % 2 == 0) else -1.0) for k in range(n)], dtype=float)


def orthonormal_residual_fro(b: np.ndarray) -> float:
    gram = b.T @ b
    return float(np.linalg.norm(gram - np.eye(gram.shape[0])))


def orthonormal_residual_max_abs(b: np.ndarray) -> float:
    gram = b.T @ b
    return float(np.max(np.abs(gram - np.eye(gram.shape[0]))))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_PROJECTIVE_CLOSURE_F672, IN_PSI_HESSIAN_F459]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P696",
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
    if not isinstance(charts, dict):
        artifact = {
            "stage": "P696",
            "status": "INVALID_F672_CHARTS_SHAPE",
            "as_of": AS_OF,
            "error": "F672 must export charts as dict",
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
            "stage": "P696",
            "status": "NOT_COMPUTABLE_MISSING_PAIR_CHARTS_IN_F672",
            "as_of": AS_OF,
            "missing_pairs": missing_pairs,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    # Load H_psi (F459)
    f459 = load_json(IN_PSI_HESSIAN_F459)
    out459 = f459.get("outputs") or {}
    h_psi = out459.get("h_psi")
    if not is_numeric_matrix(h_psi, 12):
        artifact = {
            "stage": "P696",
            "status": "INVALID_F459_H_PSI_SHAPE",
            "as_of": AS_OF,
            "error": "F459 outputs.h_psi must be a 12x12 numeric matrix",
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    H = np.array(h_psi, dtype=float)
    sym_residual_inf = float(np.max(np.abs(H - H.T)))
    Hs = 0.5 * (H + H.T)

    # Load A_m exports via F672 references.
    a_refs: dict[str, Path] = {}
    missing_a_refs: list[str] = []
    for pair in expected_pairs:
        chart = charts.get(pair)
        if not isinstance(chart, dict):
            missing_a_refs.append(f"F672.charts[{pair}] invalid")
            continue
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
            "stage": "P696",
            "status": "NOT_COMPUTABLE_MISSING_PROJECTIVE_STATE_PROJECTORS",
            "as_of": AS_OF,
            "missing": missing_a_refs,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    # Build selector-aligned channels:
    # e0, then for each pair_m: (u_m, u_m_perp), finally e6.
    n = 12
    e0 = fourier_e0(n)
    e6 = fourier_e6_on_n12()

    basis_cols: list[np.ndarray] = [e0]
    basis_labels: list[str] = ["e0"]
    pair_plane_results: dict[str, Any] = {}

    tol_u_reconstruction = 1e-10
    for pair in expected_pairs:
        m = int(pair.replace("pair", ""))
        a = load_json(a_refs[pair])
        data = a.get("data") or {}

        u_key = f"u_{m}"
        u_export = data.get(u_key)
        coords_key = f"u_{m}_coords_in_c{m}_s{m}"
        coords = data.get(coords_key)
        if not is_numeric_list_len(u_export, 12):
            artifact = {
                "stage": "P696",
                "status": "INVALID_A_M_U_VECTOR_SHAPE",
                "as_of": AS_OF,
                "error": f"{a_refs[pair].relative_to(REPO)} must contain data.{u_key} as length-12 numeric list",
                "no_false_pass": True,
            }
            OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            print(OUT_SUMMARY)
            return
        if not is_numeric_list_len(coords, 2):
            artifact = {
                "stage": "P696",
                "status": "INVALID_A_M_U_COORDS_SHAPE",
                "as_of": AS_OF,
                "error": f"{a_refs[pair].relative_to(REPO)} must contain data.{coords_key} as length-2 numeric list",
                "no_false_pass": True,
            }
            OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            print(OUT_SUMMARY)
            return

        uc = float(coords[0])
        us = float(coords[1])
        c, s = fourier_pair_basis(n, m)
        u_pred = uc * c + us * s
        u_vec = np.array([float(x) for x in u_export], dtype=float)
        u_reconstruction_max_abs_diff = float(np.max(np.abs(u_pred - u_vec)))

        if u_reconstruction_max_abs_diff > tol_u_reconstruction:
            artifact = {
                "stage": "P696",
                "status": "NOT_COMPUTABLE_FOURIER_BASIS_MISMATCH",
                "as_of": AS_OF,
                "pair": pair,
                "error": "Analytic Fourier (c_m,s_m) convention does not reconstruct exported u_m from exported coords within tolerance",
                "u_reconstruction_max_abs_diff": u_reconstruction_max_abs_diff,
                "tolerance": tol_u_reconstruction,
                "no_false_pass": True,
            }
            OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
            print(OUT_SUMMARY)
            return

        # Deterministic orthogonal complement ray in the same plane (sign-gauge tracked; ray-level use).
        v_pred = (-us) * c + (uc) * s
        u_norm = float(np.linalg.norm(u_pred))
        v_norm = float(np.linalg.norm(v_pred))
        uv_dot = float(u_pred @ v_pred)

        # Rayleigh / Ritz quantities on the symmetric H_psi.
        a_uu = float(u_pred @ (Hs @ u_pred))
        a_uv = float(u_pred @ (Hs @ v_pred))
        a_vv = float(v_pred @ (Hs @ v_pred))
        block = np.array([[a_uu, a_uv], [a_uv, a_vv]], dtype=float)
        ritz = np.linalg.eigvalsh(block)
        ritz_sorted = [float(ritz[0]), float(ritz[1])]

        # Diagonalizing angle inside the plane for the symmetric 2x2 block.
        # phi defined so that R(phi)^T block R(phi) is diagonal. (Convention: 0.5 atan2(2b, c-a))
        phi = 0.5 * math.atan2(2.0 * a_uv, (a_vv - a_uu))

        pair_plane_results[pair] = {
            "A_m_ref": str(a_refs[pair].relative_to(REPO)),
            "u_coords_in_c_s": [uc, us],
            "u_reconstruction_max_abs_diff": u_reconstruction_max_abs_diff,
            "u_norm": u_norm,
            "v_norm": v_norm,
            "u_dot_v": uv_dot,
            "block_matrix_in_(u,v)": [[a_uu, a_uv], [a_uv, a_vv]],
            "m2_proxy_u_rayleigh": a_uu,
            "m2_proxy_v_perp_rayleigh": a_vv,
            "ritz_eigenvalues_restriction_to_pair_plane": ritz_sorted,
            "diagonalizing_angle_phi": float(phi),
        }

        basis_cols.append(u_pred)
        basis_labels.append(f"{pair}:u")
        basis_cols.append(v_pred)
        basis_labels.append(f"{pair}:u_perp")

    basis_cols.append(e6)
    basis_labels.append("e6")

    B = np.column_stack(basis_cols)
    basis_orth_res_fro = orthonormal_residual_fro(B)
    basis_orth_res_max_abs = orthonormal_residual_max_abs(B)

    M = B.T @ Hs @ B
    diag = np.diag(M).copy()
    offdiag = M - np.diag(diag)
    offdiag_fro = float(np.linalg.norm(offdiag))
    diag_fro = float(np.linalg.norm(np.diag(diag)))
    offdiag_to_diag_ratio = float(offdiag_fro / max(diag_fro, 1e-300))

    # Block coupling summary for the coarse blocks: e0, pair1..pair5, e6.
    blocks = [("e0", 0, 1)]
    idx = 1
    for m in range(1, 6):
        blocks.append((f"pair{m}", idx, idx + 2))
        idx += 2
    blocks.append(("e6", 11, 12))
    offblock_max_fro = 0.0
    offblock_entries: list[dict[str, Any]] = []
    for i in range(len(blocks)):
        for j in range(i + 1, len(blocks)):
            bi, i0, i1 = blocks[i]
            bj, j0, j1 = blocks[j]
            sub = M[i0:i1, j0:j1]
            nrm = float(np.linalg.norm(sub))
            offblock_max_fro = max(offblock_max_fro, nrm)
            offblock_entries.append({"block_i": bi, "block_j": bj, "fro_norm": nrm})
    offblock_entries_sorted = sorted(offblock_entries, key=lambda e: -float(e["fro_norm"]))

    # Channel proxies (diagonal entries) as a simple operational spectrum proxy.
    channel_m2 = {lab: float(diag[i]) for i, lab in enumerate(basis_labels)}
    channel_sorted = sorted(channel_m2.items(), key=lambda kv: float(kv[1]))

    status = "PASS_SELECTOR_ALIGNED_CHANNEL_SPECTRUM_PROXY_COMPUTABLE_FROM_PROJECTIVE_SELECTOR_CLOSURE"
    artifact = {
        "stage": "P696",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": "compute_a_selector_aligned_12_channel_spectrum_proxy_for_H_psi_using_only_projective_selector_closure_rays_and_pair_plane_orthogonal_complements__no_false_pass",
        "inputs": {
            "projective_selector_closure_ref": str(IN_PROJECTIVE_CLOSURE_F672.relative_to(REPO)),
            "psi_hessian_ref": str(IN_PSI_HESSIAN_F459.relative_to(REPO)),
            "orientation_projector_operators_by_pair": {p: str(a_refs[p].relative_to(REPO)) for p in expected_pairs},
        },
        "audits": {
            "h_psi_symmetry_residual_inf_norm": sym_residual_inf,
            "u_reconstruction_tolerance_max_abs": tol_u_reconstruction,
            "selector_aligned_basis_order": list(basis_labels),
            "selector_aligned_basis_orthonormal_residual_fro": basis_orth_res_fro,
            "selector_aligned_basis_orthonormal_residual_max_abs": basis_orth_res_max_abs,
            "H_in_selector_aligned_basis_offdiag_fro": offdiag_fro,
            "H_in_selector_aligned_basis_diag_fro": diag_fro,
            "H_in_selector_aligned_basis_offdiag_to_diag_ratio": offdiag_to_diag_ratio,
            "offblock_max_fro": float(offblock_max_fro),
        },
        "results": {
            "pair_plane_restriction_results": pair_plane_results,
            "channel_m2_proxy_diagonal_of_H_in_selector_aligned_basis": channel_m2,
            "channel_m2_proxy_sorted_ascending": [{"channel": k, "m2": float(v)} for k, v in channel_sorted],
            "offblock_coupling_top10_by_fro_norm": offblock_entries_sorted[:10],
        },
        "hard_limits": [
            "Probe-only operational spectrum proxy; does not claim Standard Model identification of any computed numbers.",
            "Uses a specific diagonal/local strict-derived H_psi value instantiation (F459); does not promote it to host matching or global physics.",
            "Uses projective/ray selector data only (u and -u identified); does not claim any directed/sign-sensitive physical orientation datum in strict core.",
            "Does not claim that H_psi diagonalizes in the selector-aligned basis (mixing is audited explicitly).",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P696",
        "status": status,
        "channel_m2_proxy_sorted_ascending": [{"channel": k, "m2": float(v)} for k, v in channel_sorted],
        "H_in_selector_aligned_basis_offdiag_to_diag_ratio": offdiag_to_diag_ratio,
        "offblock_max_fro": float(offblock_max_fro),
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

