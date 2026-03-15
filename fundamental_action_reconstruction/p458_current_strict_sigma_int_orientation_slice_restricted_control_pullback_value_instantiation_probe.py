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
IN_RESTRICTED = GENERATED / "m_control_residual_diag_declared_restricted_to_sigma_int_orientation_slice_v1.json"

OUT_OBJECT = GENERATED / "m_control_residual_diag_declared_restricted_to_sigma_int_orientation_slice_value_instantiated_v1.json"
OUT_JSON = (
    GENERATED
    / "p458_current_strict_sigma_int_orientation_slice_restricted_control_pullback_value_instantiation_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p458_current_strict_sigma_int_orientation_slice_restricted_control_pullback_value_instantiation_probe_summary.json"
)


def is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def load_json_path(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def require_length_12_numeric_list(name: str, value: Any, missing: list[str]) -> list[float] | None:
    if not (isinstance(value, list) and len(value) == 12):
        missing.append(f"{name} (length-12 numeric list)")
        return None
    out: list[float] = []
    for i, x in enumerate(value):
        if not is_number(x):
            missing.append(f"{name}[{i}] (finite number)")
            out.append(float("nan"))
        else:
            out.append(float(x))
    return out


def dot(a: list[float], b: list[float]) -> float:
    return float(sum(float(x) * float(y) for x, y in zip(a, b)))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = (IN_R14, IN_R15, IN_VPSI_G4_G6, IN_RESTRICTED)
    missing_files: list[str] = [str(p.relative_to(REPO)) for p in required if not p.exists()]
    if missing_files:
        summary = {
            "stage": "P458",
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
    restricted = load_json_path(IN_RESTRICTED)

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

    vpsi = require_length_12_numeric_list("input.vpsi", inputs.get("vpsi"), missing)
    g4 = require_length_12_numeric_list("input.g4", inputs.get("g4"), missing)
    g6 = require_length_12_numeric_list("input.g6", inputs.get("g6"), missing)

    if vpsi is not None:
        for i, x in enumerate(vpsi):
            if is_number(x) and float(x) == 0.0:
                missing.append(f"input.vpsi[{i}] (must be nonzero for N477 / division by vpsi_k)")

    expected_obj = "M_control_residual_diag_declared_restricted_to_sigma_int_orientation_slice_v1"
    if restricted.get("object") != expected_obj:
        missing.append(f"P457.restricted_object (expected {expected_obj})")

    coeffs = restricted.get("restricted_coefficients_by_carrier_slot")
    if not isinstance(coeffs, dict):
        missing.append("P457.restricted_coefficients_by_carrier_slot (dict)")
        coeffs = {}

    c00 = require_length_12_numeric_list("P457.coefficients.entry_00", coeffs.get("entry_00"), missing)
    c01 = require_length_12_numeric_list("P457.coefficients.entry_01", coeffs.get("entry_01"), missing)
    c10 = require_length_12_numeric_list("P457.coefficients.entry_10", coeffs.get("entry_10"), missing)
    c11 = require_length_12_numeric_list("P457.coefficients.entry_11", coeffs.get("entry_11"), missing)

    if missing:
        summary = {
            "stage": "P458",
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
    assert c00 is not None and c01 is not None and c10 is not None and c11 is not None

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

    m00 = dot(c00, d)
    m01 = dot(c01, d)
    m10 = dot(c10, d)
    m11 = dot(c11, d)

    M = np.array([[m00, m01], [m10, m11]], dtype=float)
    eigvals = np.linalg.eigvalsh((M + M.T) * 0.5)

    # Minimal provenance hygiene: reuse the same heuristic as P437 (no stationarity enforcement).
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
        "object": "M_control_residual_diag_declared_restricted_to_sigma_int_orientation_slice_value_instantiated_v1",
        "status": status,
        "as_of": "2026-03-15",
        "intent": (
            "Value-instantiate the P457 restriction artifact (declared residual local diagonal control pullback restricted to the sigma-int orientation slice) "
            "by (i) computing a numeric diagonal residual profile d_k via the conditional N477 Yukawa-free rewrite, using the strict-frozen kernel K_total (R14), "
            "the host scalar floor m0^2 (R15), and a per-site value-provider candidate (vpsi,g4,g6), then (ii) contracting the P457 coefficient vectors with that profile. "
            "No stationarity witness is exported; the N477 rewrite remains a conditional tool."
        ),
        "inputs": {
            "r14_k_total": str(IN_R14.relative_to(REPO)),
            "r15_m0_squared": str(IN_R15.relative_to(REPO)),
            "vpsi_g4_g6_candidate": str(IN_VPSI_G4_G6.relative_to(REPO)),
            "p457_restriction_artifact": str(IN_RESTRICTED.relative_to(REPO)),
        },
        "computed": {
            "d_local_residual_profile_n477": d,
            "restricted_matrix": [[float(m00), float(m01)], [float(m10), float(m11)]],
            "symmetric_part_max_abs_diff": float(np.max(np.abs(M - M.T))),
            "trace": float(np.trace(M)),
            "det": float(np.linalg.det(M)),
            "frobenius_norm": float(np.linalg.norm(M)),
            "eigvals": [float(x) for x in eigvals.tolist()],
        },
        "hard_limits": [
            "Uses the conditional N477 Yukawa-free diagonal residual rewrite; does not export a vacuum stationarity witness.",
            "Instantiates only the P457 restriction of the declared residual local diagonal sector (D_local_residual), not the full Psi-sector Hessian block.",
            "Does not claim the instantiated restricted matrix vanishes or matches the host.",
            "Does not claim strict-core selector closure / admissible S_sel_int.",
            "Does not claim global discharge of QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    OUT_OBJECT.write_text(json.dumps(instantiated, indent=2, ensure_ascii=True) + "\n", encoding="ascii")

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P458",
        "goal": "value_instantiate_p457_sigma_int_orientation_slice_restricted_control_pullback_under_n477_conditional_diagonal_residual_rewrite",
        "inputs": instantiated["inputs"],
        "constructed_object": str(OUT_OBJECT.relative_to(REPO)),
        "provenance_hygiene": {
            "input_marked_strict_derived": bool(input_marked_strict),
            "provider_marked_strict_derived": bool(provider_marked_strict),
            "provider_theorem_refs_present": bool(provider_theorem_refs_present),
            "theorem_level_pass": bool(theorem_level_pass),
        },
        "status": status,
        "no_false_pass": True,
        "hard_limits": instantiated["hard_limits"],
    }

    summary = {
        "stage": "P458",
        "status": status,
        "constructed_object": report["constructed_object"],
        "restricted_matrix": instantiated["computed"]["restricted_matrix"],
        "frobenius_norm": instantiated["computed"]["frobenius_norm"],
        "theorem_level_pass": bool(theorem_level_pass),
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()

