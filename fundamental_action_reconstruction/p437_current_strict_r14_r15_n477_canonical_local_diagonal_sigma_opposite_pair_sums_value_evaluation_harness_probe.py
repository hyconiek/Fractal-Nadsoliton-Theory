#!/usr/bin/env python3
from __future__ import annotations

import argparse
import cmath
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

R14_JSON = GENERATED / "r14_explicit_frozen_kernel_specialization_packet_for_host_matching_route.json"
R15_JSON = GENERATED / "r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route.json"

IN_JSON = GENERATED / "p437_input_vpsi_g4_g6_candidate.json"

OUT_JSON = (
    GENERATED
    / "p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe_summary.json"
)


def is_number(x: Any) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "P437 evaluation harness for N477: compute Yukawa-free diagonal residual profile and opposite-pair sigma sums "
            "from (vpsi,g4,g6) + (R14 kernel, R15 floor)."
        )
    )
    p.add_argument("--in", dest="in_json", default=str(IN_JSON), help="Input JSON with keys: vpsi,g4,g6 (length-12 lists).")
    p.add_argument("--out-json", dest="out_json", default=str(OUT_JSON), help="Output artifact JSON path.")
    p.add_argument("--out-summary", dest="out_summary", default=str(OUT_SUMMARY), help="Output summary JSON path.")
    return p.parse_args()


def main() -> None:
    GENERATED.mkdir(exist_ok=True)
    args = parse_args()
    in_json = Path(args.in_json)
    out_json = Path(args.out_json)
    out_summary = Path(args.out_summary)

    missing: list[str] = []

    if not R14_JSON.exists():
        missing.append(str(R14_JSON.relative_to(REPO)))
    if not R15_JSON.exists():
        missing.append(str(R15_JSON.relative_to(REPO)))
    if not in_json.exists():
        missing.append(str(in_json))

    if missing:
        summary = {
            "stage": "P437",
            "status": "NOT_COMPUTABLE_MISSING_REQUIRED_FILES",
            "missing_files": missing,
            "theorem_level_pass": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }
        out_summary.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(out_summary)
        return

    r14 = load_json(R14_JSON)
    r15 = load_json(R15_JSON)
    inputs = load_json(in_json)

    K_total = r14.get("host_kernel_rows")
    if not (isinstance(K_total, list) and len(K_total) == 12 and all(isinstance(r, list) and len(r) == 12 for r in K_total)):
        missing.append("R14.host_kernel_rows (12x12 numeric K_total matrix)")

    m0_sq = (
        (r15.get("diagonal_decomposition") or {}).get("host_diagonal_floor_value")
        if isinstance(r15.get("diagonal_decomposition"), dict)
        else None
    )
    if not is_number(m0_sq):
        missing.append("R15.diagonal_decomposition.host_diagonal_floor_value (numeric m0^2)")

    vpsi = inputs.get("vpsi")
    g4 = inputs.get("g4")
    g6 = inputs.get("g6")

    def require_array(name: str, arr: Any) -> list[float] | None:
        if not (isinstance(arr, list) and len(arr) == 12):
            missing.append(f"input.{name} (length-12 numeric list)")
            return None
        out: list[float] = []
        for i, x in enumerate(arr):
            if not is_number(x):
                missing.append(f"input.{name}[{i}] (finite number)")
                out.append(float("nan"))
            else:
                out.append(float(x))
        return out

    vpsi_v = require_array("vpsi", vpsi)
    g4_v = require_array("g4", g4)
    g6_v = require_array("g6", g6)

    if vpsi_v is not None:
        for i, x in enumerate(vpsi_v):
            if is_number(x) and x == 0.0:
                missing.append(f"input.vpsi[{i}] (must be nonzero for N477 / division by vpsi_k)")

    if missing:
        summary = {
            "stage": "P437",
            "status": "NOT_COMPUTABLE_MISSING_INPUT_VALUES",
            "missing_or_invalid": sorted(set(missing)),
            "theorem_level_pass": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }
        out_summary.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(out_summary)
        return

    assert isinstance(K_total, list)
    assert isinstance(m0_sq, (int, float))
    assert vpsi_v is not None and g4_v is not None and g6_v is not None

    vpsi_abs = [abs(float(x)) for x in vpsi_v]
    min_abs_vpsi = float(min(vpsi_abs))
    max_abs_vpsi = float(max(vpsi_abs))
    vpsi_abs_ratio = float(max_abs_vpsi / min_abs_vpsi) if min_abs_vpsi > 0.0 else float("inf")

    # N477: Yukawa-free, K_total-numeric diagonal residual (conditional on stationarity + vpsi_k != 0).
    d: list[float] = []
    for k in range(12):
        denom = vpsi_v[k]
        mix = 0.0
        for j in range(12):
            if j == k:
                continue
            mix += float(K_total[k][j]) * (vpsi_v[j] / denom)
        val = -mix + 2.0 * g4_v[k] * (denom**2) + 4.0 * g6_v[k] * (denom**4) - float(m0_sq)
        d.append(float(val))

    # Opposite-pair sums (k=0..5).
    sigmas: dict[str, float] = {}
    sigma_list: list[float] = []
    for k in range(6):
        s = d[k] + d[k + 6]
        sigma_list.append(float(s))
        sigmas[f"Sigma_psi{k}_psi{k+6}"] = float(s)

    # Convenience: compute F2(d) from the six-sum reduction (N467).
    F2 = 0j
    for k in range(6):
        F2 += sigma_list[k] * cmath.exp(1j * math.pi * k / 3.0)

    artifact = {
        "stage": "P437",
        "goal": "evaluation_harness_for_N477_K_total_numeric_yukawa_free_diagonal_residual_and_opposite_pair_sums",
        "inputs": {
            "r14_k_total_source": str(R14_JSON.relative_to(REPO)),
            "r15_m0_squared_source": str(R15_JSON.relative_to(REPO)),
            "input_candidate_location": str(in_json),
        },
        "computed": {
            "vpsi_abs_min": min_abs_vpsi,
            "vpsi_abs_max": max_abs_vpsi,
            "vpsi_abs_ratio_max_over_min": vpsi_abs_ratio,
            "d_local_residual_profile": d,
            "sigma_opposite_pair_sums": sigmas,
            "F2": {"re": float(F2.real), "im": float(F2.imag), "abs": float(abs(F2))},
        },
        "provenance_note": "This is an evaluation harness only. Do not promote computed values into strict core unless inputs are strict-derived with explicit provenance.",
        "no_false_pass": True,
    }

    # Minimal provenance hygiene: if the input object is explicitly strict-derived AND it references a strict-derived
    # provider with explicit theorem anchors, then the computed six-sum instantiation is admitted as theorem-level
    # computable on the current strict export class (still not a global closure claim).
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
                prov = load_json(sp)
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

    summary = {
        "stage": "P437",
        "status": status,
        "computed_sigma_keys": list(sigmas.keys()),
        "F2_abs": float(abs(F2)),
        "vpsi_abs_min": min_abs_vpsi,
        "vpsi_abs_ratio_max_over_min": vpsi_abs_ratio,
        "input_marked_strict_derived": bool(input_marked_strict),
        "provider_marked_strict_derived": bool(provider_marked_strict),
        "provider_theorem_refs_present": bool(provider_theorem_refs_present),
        "theorem_level_pass": bool(theorem_level_pass),
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    out_json.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    out_summary.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out_json)


if __name__ == "__main__":
    main()
