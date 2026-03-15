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

IN_R14 = GENERATED / "r14_explicit_frozen_kernel_specialization_packet_for_host_matching_route.json"
IN_R15 = GENERATED / "r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route.json"
IN_VPSI_G4_G6 = GENERATED / "p437_input_vpsi_g4_g6_candidate.json"
IN_R16 = GENERATED / "r16_explicit_residual_local_diagonal_declared_control_pullback_packet_for_host_matching_route.json"

OUT_OBJECT = GENERATED / "m_control_residual_diag_declared_value_instantiated_v1.json"
OUT_JSON = GENERATED / "p459_current_strict_r16_declared_control_pullback_value_instantiation_probe.json"
OUT_SUMMARY = GENERATED / "p459_current_strict_r16_declared_control_pullback_value_instantiation_probe_summary.json"


def is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def load_json_path(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def require_length_n_numeric_list(name: str, value: Any, n: int, missing: list[str]) -> list[float] | None:
    if not (isinstance(value, list) and len(value) == n):
        missing.append(f"{name} (length-{n} numeric list)")
        return None
    out: list[float] = []
    for i, x in enumerate(value):
        if not is_number(x):
            missing.append(f"{name}[{i}] (finite number)")
            out.append(float("nan"))
        else:
            out.append(float(x))
    return out


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = (IN_R14, IN_R15, IN_VPSI_G4_G6, IN_R16)
    missing_files: list[str] = [str(p.relative_to(REPO)) for p in required if not p.exists()]
    if missing_files:
        summary = {
            "stage": "P459",
            "status": "NOT_COMPUTABLE_MISSING_REQUIRED_FILES",
            "missing_files": missing_files,
            "theorem_level_pass": False,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(json.dumps(summary, ensure_ascii=True))
        return

    r14 = load_json_path(IN_R14)
    r15 = load_json_path(IN_R15)
    inputs = load_json_path(IN_VPSI_G4_G6)
    r16 = load_json_path(IN_R16)

    missing: list[str] = []

    K_total = r14.get("host_kernel_rows")
    if not (
        isinstance(K_total, list)
        and len(K_total) == 12
        and all(isinstance(row, list) and len(row) == 12 and all(is_number(x) for x in row) for row in K_total)
    ):
        missing.append("R14.host_kernel_rows (12x12 finite numeric matrix)")

    m0_sq = (
        (r15.get("diagonal_decomposition") or {}).get("host_diagonal_floor_value")
        if isinstance(r15.get("diagonal_decomposition"), dict)
        else None
    )
    if not is_number(m0_sq):
        missing.append("R15.diagonal_decomposition.host_diagonal_floor_value (finite number)")

    vpsi = require_length_n_numeric_list("input.vpsi", inputs.get("vpsi"), 12, missing)
    g4 = require_length_n_numeric_list("input.g4", inputs.get("g4"), 12, missing)
    g6 = require_length_n_numeric_list("input.g6", inputs.get("g6"), 12, missing)
    if vpsi is not None:
        for i, x in enumerate(vpsi):
            if is_number(x) and float(x) == 0.0:
                missing.append(f"input.vpsi[{i}] (must be nonzero for N477 / division by vpsi_k)")

    control_basis = (r16.get("declared_control_pullback_packet") or {}).get("control_basis")
    if control_basis != ["c1", "s1", "c2", "s2"]:
        missing.append("R16.declared_control_pullback_packet.control_basis == ['c1','s1','c2','s2']")

    coeff_tensor = (r16.get("declared_control_pullback_packet") or {}).get("coefficient_tensor_by_carrier_slot")
    if not (isinstance(coeff_tensor, list) and len(coeff_tensor) == 4):
        missing.append("R16.declared_control_pullback_packet.coefficient_tensor_by_carrier_slot (4x4x12 numeric tensor)")
        coeff_tensor = None

    tensor: np.ndarray | None = None
    if coeff_tensor is not None:
        try:
            tensor = np.array(coeff_tensor, dtype=float)
            if tensor.shape != (4, 4, 12) or not np.isfinite(tensor).all():
                missing.append("R16.coefficient_tensor_by_carrier_slot (expected shape 4x4x12, finite)")
                tensor = None
        except Exception:
            missing.append("R16.coefficient_tensor_by_carrier_slot (numeric parse failed)")
            tensor = None

    if missing:
        summary = {
            "stage": "P459",
            "status": "NOT_COMPUTABLE_MISSING_OR_INVALID_INPUTS",
            "missing_or_invalid": sorted(set(missing)),
            "theorem_level_pass": False,
            "no_false_pass": True,
        }
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(json.dumps(summary, ensure_ascii=True))
        return

    assert isinstance(K_total, list)
    assert isinstance(m0_sq, (int, float))
    assert vpsi is not None and g4 is not None and g6 is not None
    assert tensor is not None

    # N477: Yukawa-free diagonal residual profile (conditional on stationarity + vpsi_k != 0).
    d: list[float] = []
    for k in range(12):
        denom = float(vpsi[k])
        mix = 0.0
        for j in range(12):
            if j == k:
                continue
            mix += float(K_total[k][j]) * (float(vpsi[j]) / denom)
        d_k = -mix + 2.0 * float(g4[k]) * (denom**2) + 4.0 * float(g6[k]) * (denom**4) - float(m0_sq)
        d.append(float(d_k))

    # R16: M_control_residual_diag_declared = sum_i d_i * C_i (declared residual diagonal sector only).
    d_vec = np.array(d, dtype=float)
    M = np.tensordot(tensor, d_vec, axes=([2], [0]))  # 4x4
    sym_res = float(np.max(np.abs(M - M.T)))

    pair1 = M[np.ix_([0, 1], [0, 1])]
    pair2 = M[np.ix_([2, 3], [2, 3])]
    cross_12 = M[np.ix_([0, 1], [2, 3])]
    cross_21 = M[np.ix_([2, 3], [0, 1])]

    pair1_eigs = np.linalg.eigvalsh((pair1 + pair1.T) * 0.5)
    pair2_eigs = np.linalg.eigvalsh((pair2 + pair2.T) * 0.5)

    # Minimal provenance hygiene: reuse the same heuristic as P437/P458 (no stationarity enforcement).
    input_status = str(inputs.get("status") or "")
    input_marked_strict = ("STRICT_DERIVED" in input_status.upper()) or ("strict_derived" in input_status.lower())
    provider_marked_strict = False
    provider_theorem_refs_present = False
    source_provider = inputs.get("source_provider")
    if isinstance(source_provider, str) and source_provider.strip():
        sp = Path(source_provider)
        if not sp.is_absolute():
            sp = REPO / sp
        if sp.exists():
            try:
                prov = load_json_path(sp)
                prov_status = str(prov.get("status") or "")
                prov_class = str(prov.get("classification") or "")
                provider_marked_strict = ("STRICT_DERIVED" in prov_status.upper()) or ("strict_derived" in prov_class.lower())
                theorem_refs = prov.get("theorems")
                provider_theorem_refs_present = isinstance(theorem_refs, list) and len(theorem_refs) > 0
            except Exception:
                provider_marked_strict = False
                provider_theorem_refs_present = False

    theorem_level_pass = bool(input_marked_strict and provider_marked_strict and provider_theorem_refs_present)
    status = "PASS_COMPUTED_FROM_INPUTS_REQUIRES_PROVENANCE_REVIEW"
    if theorem_level_pass:
        status = "PASS_COMPUTED_FROM_STRICT_DERIVED_INPUTS"

    instantiated = {
        "object": "M_control_residual_diag_declared_value_instantiated_v1",
        "status": status,
        "as_of": "2026-03-15",
        "intent": (
            "Value-instantiate the full declared residual local-diagonal control pullback (R16: M_control_residual_diag_declared = T_control^T D_local_residual T_control) "
            "by computing a numeric diagonal residual profile d_k via the conditional N477 Yukawa-free rewrite using (R14 K_total, R15 m0^2, and a per-site (vpsi,g4,g6) provider), "
            "then contracting the R16 coefficient tensor with that profile. No stationarity witness or host cancellation is claimed."
        ),
        "inputs": {
            "r14_k_total": str(IN_R14.relative_to(REPO)),
            "r15_m0_squared": str(IN_R15.relative_to(REPO)),
            "vpsi_g4_g6_candidate": str(IN_VPSI_G4_G6.relative_to(REPO)),
            "r16_declared_control_pullback": str(IN_R16.relative_to(REPO)),
        },
        "control_basis_order": ["c1", "s1", "c2", "s2"],
        "computed": {
            "d_local_residual_profile_n477": d,
            "m_control_residual_diag_declared": [[float(x) for x in row] for row in M.tolist()],
            "symmetric_part_max_abs_diff": sym_res,
            "pair1_block": [[float(x) for x in row] for row in pair1.tolist()],
            "pair2_block": [[float(x) for x in row] for row in pair2.tolist()],
            "cross_pair1_to_pair2_block": [[float(x) for x in row] for row in cross_12.tolist()],
            "cross_pair2_to_pair1_block": [[float(x) for x in row] for row in cross_21.tolist()],
            "pair1_eigvals": [float(x) for x in pair1_eigs.tolist()],
            "pair2_eigvals": [float(x) for x in pair2_eigs.tolist()],
            "frobenius_norm": float(np.linalg.norm(M)),
        },
        "hard_limits": [
            "Uses the conditional N477 Yukawa-free diagonal residual rewrite; does not export a vacuum stationarity witness.",
            "Instantiates only the declared residual local diagonal sector (D_local_residual), not the full Psi-sector Hessian block.",
            "Does not claim the instantiated matrix vanishes or matches the host.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim global discharge of QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    OUT_OBJECT.write_text(json.dumps(instantiated, indent=2, ensure_ascii=True) + "\n", encoding="ascii")

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P459",
        "goal": "value_instantiate_r16_declared_residual_local_diagonal_control_pullback_under_n477_conditional_diagonal_residual_rewrite",
        "inputs": instantiated["inputs"],
        "constructed_object": str(OUT_OBJECT.relative_to(REPO)),
        "provenance_hygiene": {
            "input_marked_strict_derived": bool(input_marked_strict),
            "provider_marked_strict_derived": bool(provider_marked_strict),
            "provider_theorem_refs_present": bool(provider_theorem_refs_present),
            "theorem_level_pass": bool(theorem_level_pass),
        },
        "status": status,
        "checks": {
            "control_basis_order_ok": True,
            "matrix_symmetry_max_abs_diff": sym_res,
        },
        "no_false_pass": True,
        "hard_limits": instantiated["hard_limits"],
    }

    summary = {
        "stage": "P459",
        "status": status,
        "constructed_object": report["constructed_object"],
        "pair1_block": instantiated["computed"]["pair1_block"],
        "pair2_block": instantiated["computed"]["pair2_block"],
        "cross_pair1_to_pair2_frobenius_norm": float(np.linalg.norm(cross_12)),
        "matrix_symmetry_max_abs_diff": sym_res,
        "theorem_level_pass": bool(theorem_level_pass),
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()

