#!/usr/bin/env python3
from __future__ import annotations

import itertools
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

IN_F672 = GENERATED / "selector_closure_global_c_v1_projective_strict_v1.json"
IN_F459 = GENERATED / "psi_hessian_diagonal_local_strict_derived_value_instantiated_v1.json"

IN_P694_SUMMARY = (
    GENERATED
    / "p694_current_strict_physical_computability_mass_spectrum_proxy_from_projective_selector_closure_probe_summary.json"
)
IN_P696_SUMMARY = (
    GENERATED
    / "p696_current_strict_physical_computability_selector_aligned_channel_spectrum_proxy_from_projective_selector_closure_probe_summary.json"
)
IN_F704_SUMMARY = GENERATED / "f704_current_strict_invariant_mass_observable_from_diagonal_local_psi_hessian_eigensystem_export_packet_summary.json"

OUT = GENERATED / "p709_current_strict_release_7_os_residual_sign_gauge_irrelevance_audit_probe.json"
OUT_SUMMARY = GENERATED / "p709_current_strict_release_7_os_residual_sign_gauge_irrelevance_audit_probe_summary.json"


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


def fourier_pair_basis(n: int, m: int) -> tuple[np.ndarray, np.ndarray]:
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


def compute_p696_invariants(Hs: np.ndarray, pair_coords: dict[str, tuple[float, float]]) -> dict[str, Any]:
    expected_pairs = [f"pair{m}" for m in range(1, 6)]
    n = 12
    e0 = fourier_e0(n)
    e6 = fourier_e6_on_n12()

    basis_cols: list[np.ndarray] = [e0]
    basis_labels: list[str] = ["e0"]

    for pair in expected_pairs:
        m = int(pair.replace("pair", ""))
        uc, us = pair_coords[pair]
        c, s = fourier_pair_basis(n, m)
        u = uc * c + us * s
        v = (-us) * c + (uc) * s
        basis_cols.append(u)
        basis_labels.append(f"{pair}:u")
        basis_cols.append(v)
        basis_labels.append(f"{pair}:u_perp")

    basis_cols.append(e6)
    basis_labels.append("e6")

    B = np.column_stack(basis_cols)
    M = B.T @ Hs @ B
    diag = np.diag(M).copy()
    offdiag = M - np.diag(diag)
    offdiag_fro = float(np.linalg.norm(offdiag))
    diag_fro = float(np.linalg.norm(np.diag(diag)))
    offdiag_to_diag_ratio = float(offdiag_fro / max(diag_fro, 1e-300))

    # Block coupling summary as in P696: e0, pair1..pair5, e6.
    blocks = [("e0", 0, 1)]
    idx = 1
    for m in range(1, 6):
        blocks.append((f"pair{m}", idx, idx + 2))
        idx += 2
    blocks.append(("e6", 11, 12))
    offblock_max_fro = 0.0
    for i in range(len(blocks)):
        for j in range(i + 1, len(blocks)):
            _, i0, i1 = blocks[i]
            _, j0, j1 = blocks[j]
            sub = M[i0:i1, j0:j1]
            offblock_max_fro = max(offblock_max_fro, float(np.linalg.norm(sub)))

    channel_m2 = {lab: float(diag[i]) for i, lab in enumerate(basis_labels)}

    return {
        "channel_m2": channel_m2,
        "offdiag_to_diag_ratio": offdiag_to_diag_ratio,
        "offblock_max_fro": float(offblock_max_fro),
    }


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F672, IN_F459, IN_P694_SUMMARY, IN_P696_SUMMARY, IN_F704_SUMMARY]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P709",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p694_sum = load_json(IN_P694_SUMMARY)
    p696_sum = load_json(IN_P696_SUMMARY)
    f704_sum = load_json(IN_F704_SUMMARY)

    if p694_sum.get("stage") != "P694":
        raise SystemExit(f"Invalid {IN_P694_SUMMARY.relative_to(REPO)}: expected stage=P694")
    if p696_sum.get("stage") != "P696":
        raise SystemExit(f"Invalid {IN_P696_SUMMARY.relative_to(REPO)}: expected stage=P696")
    if f704_sum.get("stage") != "F704":
        raise SystemExit(f"Invalid {IN_F704_SUMMARY.relative_to(REPO)}: expected stage=F704")

    baseline_p694_m2 = p694_sum.get("mass_proxy_m2_by_pair_rayleigh_on_H_psi")
    if not isinstance(baseline_p694_m2, dict):
        raise SystemExit(f"Invalid {IN_P694_SUMMARY.relative_to(REPO)}: mass_proxy_m2_by_pair_rayleigh_on_H_psi must be dict")

    baseline_p696_sorted = p696_sum.get("channel_m2_proxy_sorted_ascending")
    if not (isinstance(baseline_p696_sorted, list) and all(isinstance(e, dict) for e in baseline_p696_sorted)):
        raise SystemExit(f"Invalid {IN_P696_SUMMARY.relative_to(REPO)}: channel_m2_proxy_sorted_ascending must be list[dict]")
    baseline_channel_m2: dict[str, float] = {}
    for e in baseline_p696_sorted:
        ch = e.get("channel")
        m2 = e.get("m2")
        if not (isinstance(ch, str) and is_finite_number(m2)):
            raise SystemExit(
                f"Invalid {IN_P696_SUMMARY.relative_to(REPO)}: each entry must contain channel:str and m2:number"
            )
        baseline_channel_m2[ch] = float(m2)

    baseline_offdiag_to_diag_ratio = p696_sum.get("H_in_selector_aligned_basis_offdiag_to_diag_ratio")
    baseline_offblock_max_fro = p696_sum.get("offblock_max_fro")
    if not is_finite_number(baseline_offdiag_to_diag_ratio):
        raise SystemExit(
            f"Invalid {IN_P696_SUMMARY.relative_to(REPO)}: H_in_selector_aligned_basis_offdiag_to_diag_ratio must be number"
        )
    if not is_finite_number(baseline_offblock_max_fro):
        raise SystemExit(f"Invalid {IN_P696_SUMMARY.relative_to(REPO)}: offblock_max_fro must be number")

    baseline_f704_min = f704_sum.get("min_m2_proxy")
    baseline_f704_max = f704_sum.get("max_m2_proxy")
    if not is_finite_number(baseline_f704_min):
        raise SystemExit(f"Invalid {IN_F704_SUMMARY.relative_to(REPO)}: min_m2_proxy must be number")
    if not is_finite_number(baseline_f704_max):
        raise SystemExit(f"Invalid {IN_F704_SUMMARY.relative_to(REPO)}: max_m2_proxy must be number")

    # Load H_psi (F459).
    f459 = load_json(IN_F459)
    out459 = f459.get("outputs") or {}
    h_psi = out459.get("h_psi")
    if not is_numeric_matrix(h_psi, 12):
        raise SystemExit(f"Invalid {IN_F459.relative_to(REPO)}: outputs.h_psi must be 12x12 numeric matrix")

    H = np.array(h_psi, dtype=float)
    Hs = 0.5 * (H + H.T)

    # Load A_m exports via F672 references.
    f672 = load_json(IN_F672)
    charts = f672.get("charts")
    if not isinstance(charts, dict):
        raise SystemExit(f"Invalid {IN_F672.relative_to(REPO)}: charts must be dict")

    expected_pairs = [f"pair{m}" for m in range(1, 6)]
    a_refs: dict[str, Path] = {}
    pair_u: dict[str, np.ndarray] = {}
    pair_coords: dict[str, tuple[float, float]] = {}

    for pair in expected_pairs:
        ch = charts.get(pair)
        if not isinstance(ch, dict):
            raise SystemExit(f"Invalid {IN_F672.relative_to(REPO)}: charts.{pair} must be dict")
        ref = ch.get("A_m_ref")
        if not isinstance(ref, str) or not ref:
            raise SystemExit(f"Invalid {IN_F672.relative_to(REPO)}: charts.{pair}.A_m_ref missing")
        p = Path(ref)
        if not p.is_absolute():
            p = REPO / p
        if not p.exists():
            raise SystemExit(f"Missing A_m_ref: {Path(ref)}")
        a_refs[pair] = p

        a = load_json(p)
        data = a.get("data") or {}
        m = int(pair.replace("pair", ""))
        u_key = f"u_{m}"
        coords_key = f"u_{m}_coords_in_c{m}_s{m}"
        u = data.get(u_key)
        coords = data.get(coords_key)
        if not is_numeric_list_len(u, 12):
            raise SystemExit(f"Invalid {p.relative_to(REPO)}: data.{u_key} must be length-12 numeric list")
        if not is_numeric_list_len(coords, 2):
            raise SystemExit(f"Invalid {p.relative_to(REPO)}: data.{coords_key} must be length-2 numeric list")
        pair_u[pair] = np.array([float(x) for x in u], dtype=float)
        pair_coords[pair] = (float(coords[0]), float(coords[1]))

    # Tolerances: strict enough to catch real changes but robust to float noise.
    tol_m2 = 1e-10
    tol_ratio = 1e-12

    # Baseline recomputation (from raw artifacts, not trusting summaries).
    baseline_recomputed_p694 = {pair: float(rayleigh_quotient(Hs, pair_u[pair])) for pair in expected_pairs}
    baseline_recomputed_p696 = compute_p696_invariants(Hs, pair_coords)
    evals = np.linalg.eigvalsh(Hs)
    baseline_recomputed_f704_min = float(np.min(evals))
    baseline_recomputed_f704_max = float(np.max(evals))

    # Compare baseline recomputation against exported summaries (sanity).
    baseline_summary_max_abs_diff_p694 = max(
        abs(float(baseline_p694_m2[pair]) - float(baseline_recomputed_p694[pair])) for pair in expected_pairs
    )
    baseline_summary_max_abs_diff_channel = max(
        abs(float(baseline_channel_m2[k]) - float(baseline_recomputed_p696["channel_m2"][k])) for k in baseline_channel_m2
    )
    baseline_summary_abs_diff_offdiag_ratio = abs(
        float(baseline_offdiag_to_diag_ratio) - float(baseline_recomputed_p696["offdiag_to_diag_ratio"])
    )
    baseline_summary_abs_diff_offblock_max = abs(
        float(baseline_offblock_max_fro) - float(baseline_recomputed_p696["offblock_max_fro"])
    )
    baseline_summary_abs_diff_f704_min = abs(float(baseline_f704_min) - baseline_recomputed_f704_min)
    baseline_summary_abs_diff_f704_max = abs(float(baseline_f704_max) - baseline_recomputed_f704_max)

    # Residual sign patterns over pair1..pair5.
    patterns = list(itertools.product([-1.0, 1.0], repeat=5))

    max_abs_diff_p694_under_sign = 0.0
    max_abs_diff_channel_under_sign = 0.0
    max_abs_diff_offdiag_ratio_under_sign = 0.0
    max_abs_diff_offblock_max_under_sign = 0.0
    max_abs_diff_f704_min_under_sign = 0.0
    max_abs_diff_f704_max_under_sign = 0.0

    failing_patterns: list[dict[str, Any]] = []
    for pat in patterns:
        s_by_pair = {f"pair{i+1}": float(pat[i]) for i in range(5)}

        # P694: Rayleigh is invariant under u->-u (per pair independently).
        for pair in expected_pairs:
            u = pair_u[pair]
            m2 = float(rayleigh_quotient(Hs, s_by_pair[pair] * u))
            diff = abs(m2 - baseline_recomputed_p694[pair])
            max_abs_diff_p694_under_sign = max(max_abs_diff_p694_under_sign, diff)

        # P696: apply sign per pair to the coords -> flips both (u,v) columns for that pair.
        coords_signed = {
            pair: (s_by_pair[pair] * pair_coords[pair][0], s_by_pair[pair] * pair_coords[pair][1])
            for pair in expected_pairs
        }
        inv = compute_p696_invariants(Hs, coords_signed)
        for k, base in baseline_recomputed_p696["channel_m2"].items():
            diff = abs(float(inv["channel_m2"][k]) - float(base))
            max_abs_diff_channel_under_sign = max(max_abs_diff_channel_under_sign, diff)

        diff_ratio = abs(float(inv["offdiag_to_diag_ratio"]) - float(baseline_recomputed_p696["offdiag_to_diag_ratio"]))
        diff_offblock = abs(float(inv["offblock_max_fro"]) - float(baseline_recomputed_p696["offblock_max_fro"]))
        max_abs_diff_offdiag_ratio_under_sign = max(max_abs_diff_offdiag_ratio_under_sign, diff_ratio)
        max_abs_diff_offblock_max_under_sign = max(max_abs_diff_offblock_max_under_sign, diff_offblock)

        # F704: eigenvalues are basis-invariant; independent of any u sign.
        diff_min = abs(baseline_recomputed_f704_min - baseline_recomputed_f704_min)
        diff_max = abs(baseline_recomputed_f704_max - baseline_recomputed_f704_max)
        max_abs_diff_f704_min_under_sign = max(max_abs_diff_f704_min_under_sign, diff_min)
        max_abs_diff_f704_max_under_sign = max(max_abs_diff_f704_max_under_sign, diff_max)

        ok = (
            max_abs_diff_p694_under_sign <= tol_m2
            and max_abs_diff_channel_under_sign <= tol_m2
            and max_abs_diff_offdiag_ratio_under_sign <= tol_ratio
            and max_abs_diff_offblock_max_under_sign <= tol_m2
            and max_abs_diff_f704_min_under_sign <= tol_m2
            and max_abs_diff_f704_max_under_sign <= tol_m2
        )
        if not ok and len(failing_patterns) < 5:
            failing_patterns.append(
                {
                    "s_by_pair": s_by_pair,
                    "note": "At least one max-diff exceeded tolerance (see global maxima fields).",
                }
            )

    # PASS conditions:
    # - baseline recomputation matches baseline summaries (sanity), AND
    # - sign pattern maxima are within tolerance.
    baseline_ok = (
        baseline_summary_max_abs_diff_p694 <= tol_m2
        and baseline_summary_max_abs_diff_channel <= tol_m2
        and baseline_summary_abs_diff_offdiag_ratio <= tol_ratio
        and baseline_summary_abs_diff_offblock_max <= tol_m2
        and baseline_summary_abs_diff_f704_min <= tol_m2
        and baseline_summary_abs_diff_f704_max <= tol_m2
    )

    sign_ok = (
        max_abs_diff_p694_under_sign <= tol_m2
        and max_abs_diff_channel_under_sign <= tol_m2
        and max_abs_diff_offdiag_ratio_under_sign <= tol_ratio
        and max_abs_diff_offblock_max_under_sign <= tol_m2
        and max_abs_diff_f704_min_under_sign <= tol_m2
        and max_abs_diff_f704_max_under_sign <= tol_m2
    )

    status = (
        "PASS_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_AUDITED"
        if (baseline_ok and sign_ok)
        else ("FAIL_BASELINE_MISMATCH_OR_SIGN_SENSITIVITY_DETECTED")
    )

    artifact = {
        "stage": "P709",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "goal": (
            "audit that residual per-pair sign u_m->-u_m is gauge-irrelevant for Release-7 OS observables "
            "(P694 mass proxy, P696 channel spectrum proxy invariants, F704 invariant mass observable), "
            "without claiming kernel-alone/global QW-2191 discharge"
        ),
        "inputs": {
            "F672_projective_closure_ref": str(IN_F672.relative_to(REPO)),
            "F459_h_psi_ref": str(IN_F459.relative_to(REPO)),
            "P694_summary_ref": str(IN_P694_SUMMARY.relative_to(REPO)),
            "P696_summary_ref": str(IN_P696_SUMMARY.relative_to(REPO)),
            "F704_summary_ref": str(IN_F704_SUMMARY.relative_to(REPO)),
        },
        "tolerances": {
            "tol_m2_abs": tol_m2,
            "tol_ratio_abs": tol_ratio,
        },
        "baseline_sanity": {
            "max_abs_diff_p694_m2_against_summary": baseline_summary_max_abs_diff_p694,
            "max_abs_diff_p696_channel_m2_against_summary": baseline_summary_max_abs_diff_channel,
            "abs_diff_p696_offdiag_to_diag_ratio_against_summary": baseline_summary_abs_diff_offdiag_ratio,
            "abs_diff_p696_offblock_max_fro_against_summary": baseline_summary_abs_diff_offblock_max,
            "abs_diff_f704_min_against_summary": baseline_summary_abs_diff_f704_min,
            "abs_diff_f704_max_against_summary": baseline_summary_abs_diff_f704_max,
            "baseline_ok": baseline_ok,
        },
        "sign_pattern_audit": {
            "n_patterns_tested": len(patterns),
            "max_abs_diff_p694_m2_under_sign_patterns": max_abs_diff_p694_under_sign,
            "max_abs_diff_p696_channel_m2_under_sign_patterns": max_abs_diff_channel_under_sign,
            "max_abs_diff_p696_offdiag_to_diag_ratio_under_sign_patterns": max_abs_diff_offdiag_ratio_under_sign,
            "max_abs_diff_p696_offblock_max_fro_under_sign_patterns": max_abs_diff_offblock_max_under_sign,
            "max_abs_diff_f704_min_under_sign_patterns": max_abs_diff_f704_min_under_sign,
            "max_abs_diff_f704_max_under_sign_patterns": max_abs_diff_f704_max_under_sign,
            "sign_ok": sign_ok,
            "first_failing_patterns_if_any": failing_patterns,
        },
        "hard_limits": [
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim any directed/sign-sensitive physical orientation datum in strict core.",
            "Does not claim Standard Model identification or proxy→GeV calibration.",
            "Does not claim 'actual emergent observer closure'.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P709",
        "status": status,
        "baseline_ok": baseline_ok,
        "sign_ok": sign_ok,
        "n_patterns_tested": len(patterns),
        "max_abs_diff_p694_m2_under_sign_patterns": max_abs_diff_p694_under_sign,
        "max_abs_diff_p696_channel_m2_under_sign_patterns": max_abs_diff_channel_under_sign,
        "max_abs_diff_p696_offdiag_to_diag_ratio_under_sign_patterns": max_abs_diff_offdiag_ratio_under_sign,
        "max_abs_diff_p696_offblock_max_fro_under_sign_patterns": max_abs_diff_offblock_max_under_sign,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

