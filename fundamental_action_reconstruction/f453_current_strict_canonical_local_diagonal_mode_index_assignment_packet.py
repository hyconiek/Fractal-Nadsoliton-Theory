#!/usr/bin/env python3
from __future__ import annotations

import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

P437_SCRIPT = ROOT / "p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe.py"
P449_SCRIPT = ROOT / "p449_current_strict_canonical_local_diagonal_multi_pair_o2_cut_defect_evaluation_probe.py"

P437_OUT = (
    GENERATED
    / "p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe.json"
)
P437_SUMMARY = (
    GENERATED
    / "p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe_summary.json"
)
P449_OUT = GENERATED / "p449_current_strict_canonical_local_diagonal_multi_pair_o2_cut_defect_evaluation_probe.json"
P449_SUMMARY = GENERATED / "p449_current_strict_canonical_local_diagonal_multi_pair_o2_cut_defect_evaluation_probe_summary.json"

QW2190 = REPO / "report_qw2190_kernel_mode_representation_emergence_gate.json"

OUT_ASSIGNMENT = GENERATED / "mode_index_assignment_canonical_local_diagonal_strict_derived_v1.json"
OUT_ASSIGNMENT_SUMMARY = GENERATED / "mode_index_assignment_canonical_local_diagonal_strict_derived_v1_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def run_script(path: Path) -> dict[str, Any]:
    proc = subprocess.run(
        ["python3", str(path)],
        cwd=str(REPO),
        capture_output=True,
        text=True,
        check=False,
    )
    return {
        "script": str(path.relative_to(REPO)),
        "returncode": proc.returncode,
        "stdout_tail": "\n".join(proc.stdout.splitlines()[-10:]),
        "stderr_tail": "\n".join(proc.stderr.splitlines()[-10:]),
    }


def real_fourier_basis(n: int) -> dict[str, np.ndarray]:
    x = np.arange(n, dtype=float)
    basis: dict[str, np.ndarray] = {"e0": np.ones(n, dtype=float) / math.sqrt(n)}
    for m in range(1, n // 2):
        basis[f"c{m}"] = math.sqrt(2.0 / n) * np.cos(2.0 * math.pi * m * x / n)
        basis[f"s{m}"] = math.sqrt(2.0 / n) * np.sin(2.0 * math.pi * m * x / n)
    if n % 2 == 0:
        basis[f"e{n//2}"] = ((-1.0) ** x) / math.sqrt(n)
    return basis


def orthonormal_residual(b: np.ndarray) -> float:
    m = b.T @ b
    return float(np.linalg.norm(m - np.eye(m.shape[0])))


def pair_basis_from_theta(c: np.ndarray, s: np.ndarray, theta: float) -> tuple[np.ndarray, np.ndarray]:
    ct = math.cos(theta)
    st = math.sin(theta)
    u_plus = ct * c + st * s
    u_minus = -st * c + ct * s
    return u_plus, u_minus


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    runs = [run_script(P437_SCRIPT), run_script(P449_SCRIPT)]
    if any(int(r["returncode"]) != 0 for r in runs):
        raise SystemExit(json.dumps({"status": "DEPENDENCY_SCRIPT_FAILED", "runs": runs}, ensure_ascii=True))

    p437_summary = load_json(P437_SUMMARY)
    p449 = load_json(P449_OUT)
    p449_summary = load_json(P449_SUMMARY)
    q2190 = load_json(QW2190)

    n = int(p449.get("computed", {}).get("n") or q2190["mode_mapping"]["n_octaves"])
    fourier = real_fourier_basis(n)
    e0 = fourier["e0"]
    e_half = fourier.get(f"e{n//2}") if (n % 2 == 0) else None

    strict_derived_inputs_ok = bool(p437_summary.get("input_marked_strict_derived") is True) and bool(
        p437_summary.get("provider_marked_strict_derived") is True
    )
    strict_derived_provenance_ok = bool(p437_summary.get("provider_theorem_refs_present") is True) and bool(
        p449_summary.get("input_theorem_level_pass_from_P437") is True
    )
    theorem_level_pass = bool(p449_summary.get("theorem_level_pass"))

    tol = 1e-12
    pairs: dict[str, Any] = {}
    basis_vectors: list[np.ndarray] = [e0]

    all_pairs_cut = True
    for m in range(1, n // 2):
        pair_key = f"pair{m}"
        pair = (p449.get("computed") or {}).get("pairs", {}).get(pair_key)
        if not isinstance(pair, dict):
            raise SystemExit(
                json.dumps(
                    {"status": "MISSING_PAIR_DEFECT", "expected": pair_key, "available": list(((p449.get("computed") or {}).get("pairs") or {}).keys())},
                    ensure_ascii=True,
                )
            )
        theta_star = float(pair["theta_star"])
        F = pair["F"]
        absF = float(F["abs"])
        cuts = absF > tol
        all_pairs_cut = all_pairs_cut and cuts

        c = fourier[f"c{m}"]
        s = fourier[f"s{m}"]
        u_plus, u_minus = pair_basis_from_theta(c, s, theta_star)
        b = np.column_stack([u_plus, u_minus])

        pairs[pair_key] = {
            "m": int(pair["m"]),
            "l": int(pair["l"]),
            "F_2m": F,
            "theta_star": theta_star,
            "cuts_O2_to_Z2_by_absF_nonzero": bool(cuts),
            "eigenvalues": pair.get("eigenvalues"),
            "u_plus": [float(x) for x in u_plus.tolist()],
            "u_minus": [float(x) for x in u_minus.tolist()],
            "audits": {
                "pair_basis_orthonormal_residual": float(np.linalg.norm(b.T @ b - np.eye(2))),
                "pair_basis_det_gram": float(np.linalg.det(b.T @ b)),
                "pair_basis_rank": int(np.linalg.matrix_rank(b, tol=1e-12)),
                "pair_basis_orthonormal_residual_vs_identity": orthonormal_residual(b),
            },
        }
        basis_vectors.extend([u_plus, u_minus])

    if e_half is not None:
        basis_vectors.append(e_half)

    full_basis = np.column_stack(basis_vectors)
    full_basis_orth_res = orthonormal_residual(full_basis)

    assignment = {
        "object": "ModeIndexAssignment_canonical_local_diagonal_strict_derived_v1",
        "status": "actual_exported_strict_core_mode_index_assignment__strict_derived_value_instantiation__no_false_pass",
        "as_of": "2026-03-15",
        "intent": (
            "Export one strict-derived canonical mode-index assignment basis on the QW-2190 real Fourier scaffold for n=12, "
            "constructed by diagonal/local defect reconstruction on each Fourier-degenerate pair plane pair_m (m=1..5) from the "
            "canonical diagonal/local residual profile computed by P437 and the multi-pair Fourier defect evaluation of P449. "
            "This provides a canonical eigenbasis representative (up to residual Z2 sign) on each degenerate pair plane in the "
            "declared diagonal/local lane (N484/N485/N487), without implying global selector closure or ToE closure."
        ),
        "scope": "strict_core_strict_derived_value_instantiation_only",
        "derived_from": {
            "p437_profile_artifact": str(P437_OUT.relative_to(REPO)),
            "p437_summary": str(P437_SUMMARY.relative_to(REPO)),
            "p449_defect_artifact": str(P449_OUT.relative_to(REPO)),
            "p449_summary": str(P449_SUMMARY.relative_to(REPO)),
            "qw2190_mode_scaffold": str(QW2190.relative_to(REPO)),
        },
        "provenance_checks": {
            "strict_derived_inputs_marked": strict_derived_inputs_ok,
            "strict_derived_provenance_chain_present": strict_derived_provenance_ok,
            "input_theorem_level_pass": theorem_level_pass,
            "note": "This packet trusts only the strict-derived diagonal/local lane artifacts (P437/P449).",
        },
        "construction": {
            "diagonal_sector_pair_m_reconstruction": {
                "theta_star_rule": "theta_* = (1/2) atan2(Im(F_2m), Re(F_2m))",
                "u_plus_rule": "u_+ := cos(theta_*) c_m + sin(theta_*) s_m",
                "u_minus_rule": "u_- := -sin(theta_*) c_m + cos(theta_*) s_m",
                "references": ["N484", "N485", "N487"],
            }
        },
        "outputs": {
            "n": n,
            "basis_vectors_order": (
                ["e0", *[f"pair{m}:u_plus,u_minus" for m in range(1, n // 2)], (f"e{n//2}" if (e_half is not None) else None)]
            ),
            "e0": [float(x) for x in e0.tolist()],
            (f"e{n//2}" if e_half is not None else "e_half_absent"): (
                [float(x) for x in e_half.tolist()] if e_half is not None else None
            ),
            "pairs": pairs,
            "all_pairs_cut": bool(all_pairs_cut),
        },
        "audits": {
            "full_basis_shape": [int(full_basis.shape[0]), int(full_basis.shape[1])],
            "full_basis_orthonormal_residual_vs_identity": full_basis_orth_res,
            "full_basis_det_gram": float(np.linalg.det(full_basis.T @ full_basis)),
            "dependency_runs": runs,
        },
        "hard_limits": [
            "Does not claim axiom-free global physical uniqueness beyond the declared diagonal/local lane and n=12 scope.",
            "Does not claim strict-core selector closure nor admissible S_sel_int.",
            "Does not claim global QW-2191 discharge (kernel-alone obstruction remains true).",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F453",
        "status": "F453_EXECUTED_CURRENT_STRICT_CANONICAL_LOCAL_DIAGONAL_MODE_INDEX_ASSIGNMENT_PACKET_NO_FALSE_PASS",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": {
            "n": n,
            "all_pairs_cut": bool(all_pairs_cut),
            "full_basis_orthonormal_residual_vs_identity": full_basis_orth_res,
        },
        "inputs": {
            "p437_summary": str(P437_SUMMARY.relative_to(REPO)),
            "p449_summary": str(P449_SUMMARY.relative_to(REPO)),
        },
        "no_false_pass": True,
    }

    OUT_ASSIGNMENT.write_text(json.dumps(assignment, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_ASSIGNMENT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_ASSIGNMENT_SUMMARY)


if __name__ == "__main__":
    main()

