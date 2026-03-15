#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F459 = GENERATED / "psi_hessian_diagonal_local_strict_derived_value_instantiated_v1.json"
IN_F453 = GENERATED / "mode_index_assignment_canonical_local_diagonal_strict_derived_v1.json"
IN_F454 = GENERATED / "mode_index_assignment_shannon_element_order_reference_strict_core_v1.json"

OUT_JSON = (
    GENERATED
    / "p463_current_strict_diagonal_local_psi_hessian_eigenmodes_to_mode_index_assignment_projection_audit_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p463_current_strict_diagonal_local_psi_hessian_eigenmodes_to_mode_index_assignment_projection_audit_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_full_basis_from_assignment(assignment: dict[str, Any]) -> tuple[np.ndarray, list[str]]:
    out = (assignment.get("outputs") or {})
    order = out.get("basis_vectors_order")
    if not isinstance(order, list):
        raise ValueError("assignment.outputs.basis_vectors_order missing or invalid")

    cols: list[np.ndarray] = []
    labels: list[str] = []

    for item in order:
        if not isinstance(item, str):
            raise ValueError("basis_vectors_order must be list[str]")
        if item == "e0":
            v = out.get("e0")
            if not isinstance(v, list) or len(v) != 12:
                raise ValueError("assignment.outputs.e0 must be length-12 list")
            cols.append(np.array([float(x) for x in v], dtype=float))
            labels.append("e0")
            continue
        if item == "e6":
            v = out.get("e6")
            if not isinstance(v, list) or len(v) != 12:
                raise ValueError("assignment.outputs.e6 must be length-12 list")
            cols.append(np.array([float(x) for x in v], dtype=float))
            labels.append("e6")
            continue

        # pair entry like "pairm:u_plus,u_minus"
        if ":" not in item:
            raise ValueError(f"unexpected basis order entry: {item}")
        pair, rest = item.split(":", 1)
        if "," not in rest:
            raise ValueError(f"unexpected pair basis spec: {item}")
        # we ignore tokens and fetch stored vectors directly
        pairs = (out.get("pairs") or {})
        pdata = pairs.get(pair)
        if not isinstance(pdata, dict):
            raise ValueError(f"assignment.outputs.pairs[{pair}] missing")
        u_plus = pdata.get("u_plus")
        u_minus = pdata.get("u_minus")
        if not (isinstance(u_plus, list) and isinstance(u_minus, list) and len(u_plus) == 12 and len(u_minus) == 12):
            raise ValueError(f"assignment.outputs.pairs[{pair}].u_plus/u_minus must be length-12 lists")
        cols.append(np.array([float(x) for x in u_plus], dtype=float))
        labels.append(f"{pair}:u_plus")
        cols.append(np.array([float(x) for x in u_minus], dtype=float))
        labels.append(f"{pair}:u_minus")

    if len(cols) != 12:
        raise ValueError(f"expected 12 basis vectors after expansion, got {len(cols)}")
    basis = np.column_stack(cols)
    return basis, labels


def orthonormal_residual(b: np.ndarray) -> float:
    gram = b.T @ b
    return float(np.linalg.norm(gram - np.eye(gram.shape[0])))


def top_k_components(weights: np.ndarray, labels: list[str], k: int = 5) -> list[dict[str, Any]]:
    idx = np.argsort(-weights)[:k]
    out = []
    for i in idx:
        out.append({"label": labels[int(i)], "weight": float(weights[int(i)])})
    return out


def pair_mass_from_weights(weights: np.ndarray, labels: list[str]) -> dict[str, float]:
    masses: dict[str, float] = {"e0": 0.0, "e6": 0.0}
    for m in range(1, 6):
        masses[f"pair{m}"] = 0.0
    for w, lab in zip(weights.tolist(), labels):
        if lab in ("e0", "e6"):
            masses[lab] += float(w)
        else:
            pair = lab.split(":", 1)[0]
            masses[pair] = masses.get(pair, 0.0) + float(w)
    return masses


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    missing: list[str] = []
    for p in (IN_F459, IN_F453, IN_F454):
        if not p.exists():
            missing.append(str(p.relative_to(REPO)))
    if missing:
        summary = {
            "stage": "P463",
            "status": "NOT_COMPUTABLE_MISSING_INPUTS",
            "missing": missing,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f459 = load_json(IN_F459)
    f453 = load_json(IN_F453)
    f454 = load_json(IN_F454)

    evals = (f459.get("outputs") or {}).get("eigenvalues_ascending")
    ecols = (f459.get("outputs") or {}).get("eigenvectors_columns")
    if not (isinstance(evals, list) and isinstance(ecols, list) and len(evals) == 12 and len(ecols) == 12):
        raise SystemExit(json.dumps({"status": "BAD_F459_EIGENSYSTEM_SHAPE", "path": str(IN_F459.relative_to(REPO))}))
    E = np.column_stack([np.array([float(x) for x in col], dtype=float) for col in ecols])
    evals_v = np.array([float(x) for x in evals], dtype=float)

    B_diag, labels_diag = build_full_basis_from_assignment(f453)
    B_sha, labels_sha = build_full_basis_from_assignment(f454)

    audits = {
        "E_orthonormal_residual_l2": orthonormal_residual(E),
        "B_diag_orthonormal_residual_l2": orthonormal_residual(B_diag),
        "B_sha_orthonormal_residual_l2": orthonormal_residual(B_sha),
    }

    O_diag = B_diag.T @ E
    O_sha = B_sha.T @ E

    abs2_diag = np.abs(O_diag) ** 2
    abs2_sha = np.abs(O_sha) ** 2

    # Per-eigenmode summaries
    per_mode: list[dict[str, Any]] = []
    for j in range(12):
        wdiag = abs2_diag[:, j]
        wsha = abs2_sha[:, j]
        entry = {
            "j": int(j),
            "lambda": float(evals_v[j]),
            "diag_basis": {
                "top_components": top_k_components(wdiag, labels_diag, k=6),
                "max_weight": float(np.max(wdiag)),
                "pair_mass": pair_mass_from_weights(wdiag, labels_diag),
                "effective_support_exp_entropy": float(math.exp(-float(np.sum(wdiag * np.log(np.maximum(wdiag, 1e-300)))))),
            },
            "shannon_basis": {
                "top_components": top_k_components(wsha, labels_sha, k=6),
                "max_weight": float(np.max(wsha)),
                "pair_mass": pair_mass_from_weights(wsha, labels_sha),
                "effective_support_exp_entropy": float(math.exp(-float(np.sum(wsha * np.log(np.maximum(wsha, 1e-300)))))),
            },
        }
        per_mode.append(entry)

    artifact = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P463",
        "goal": "audit_projection_of_f459_psi_hessian_eigenmodes_into_exported_mode_index_assignment_bases_f453_f454",
        "inputs": {
            "f459_eigensystem": str(IN_F459.relative_to(REPO)),
            "f453_diag_local_assignment": str(IN_F453.relative_to(REPO)),
            "f454_shannon_assignment": str(IN_F454.relative_to(REPO)),
        },
        "audits": audits,
        "results": {
            "overlap_abs2_diag_basis": [[float(x) for x in row.tolist()] for row in abs2_diag],
            "overlap_abs2_shannon_basis": [[float(x) for x in row.tolist()] for row in abs2_sha],
            "per_eigenmode": per_mode,
        },
        "hard_limits": [
            "Probe-only: does not claim either mode-index assignment basis diagonalizes the full H_psi operator.",
            "Does not export any global selector transition/gluing object.",
            "Does not claim any global discharge of QW-2191 beyond the declared lane-scoped instantiations.",
            "Does not claim strict-core selector closure nor admissible S_sel_int.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    # Compact summary: for each eigenmode, record dominant pair mass + max basis weight
    summary_rows = []
    for m in per_mode:
        pm = m["diag_basis"]["pair_mass"]
        top_pair = max(pm.items(), key=lambda kv: kv[1])
        summary_rows.append(
            {
                "j": m["j"],
                "lambda": m["lambda"],
                "diag_top_pair": {"label": top_pair[0], "mass": float(top_pair[1])},
                "diag_max_basis_weight": float(m["diag_basis"]["max_weight"]),
            }
        )

    summary = {
        "stage": "P463",
        "status": "PASS_PROBE_READY",
        "audits": audits,
        "per_eigenmode_compact": summary_rows,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

