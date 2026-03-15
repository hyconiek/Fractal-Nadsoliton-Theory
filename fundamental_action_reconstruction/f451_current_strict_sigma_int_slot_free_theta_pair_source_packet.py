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

IN_ALPHA_GEO = GENERATED / "alpha_geo_strict_derived_v1.json"
IN_SIGMA_INT = GENERATED / "sigma_int_strict_derived_v1.json"
IN_R_ORD = GENERATED / "r_ord_z12_v1_reference_distribution.json"
IN_QW2190 = REPO / "report_qw2190_kernel_mode_representation_emergence_gate.json"

OUT_THETA_PAIR = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"
OUT_THETA_PAIR_SUMMARY = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def is_int_sign(x: Any) -> bool:
    return isinstance(x, int) and x in (-1, 1)


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


def zn_element_order(n: int, k: int) -> int:
    kk = int(k) % int(n)
    if kk == 0:
        return 1
    return int(n) // math.gcd(kk, int(n))


def ensure_r_ord_present(alpha_geo_symbolic: str) -> dict[str, Any]:
    if IN_R_ORD.exists():
        return {"created": False, "path": str(IN_R_ORD.relative_to(REPO))}

    n = 12
    ord_by_x = [zn_element_order(n, x) for x in range(n)]
    payload = {
        "object_id": "r_ord_z12_v1_reference_distribution",
        "type": "z12_reference_distribution_v1",
        "definition": "r_ord(x) ∝ exp(-alpha_geo * ord_Z12(x))",
        "alpha_geo": {"source_object_id": "alpha_geo_strict_derived_v1", "value_symbolic": alpha_geo_symbolic},
        "carrier": {"group_object_id": "z_12_v1_group", "index_set_object_id": "i_12_v1_index_set"},
        "ord_z12_by_x": ord_by_x,
        "normalization": {
            "symbolic": "Z = Σ_{x∈Z_12} exp(-alpha_geo*ord_z12[x])",
            "notes": ["This artifact intentionally stores ord_z12 values; numeric normalization can be computed downstream."],
        },
        "invariance_notes": [
            "Aut(Z_12)-invariant reference shape (N479): no marked generator/direction.",
            "Not translation-invariant on the regular action: it distinguishes the identity orbit {0}.",
        ],
        "provenance": {"packet": "F451", "theorems": ["N479", "N480", "N488"]},
    }
    IN_R_ORD.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    return {"created": True, "path": str(IN_R_ORD.relative_to(REPO))}


def opposite_pair_sums(values: list[int]) -> list[int]:
    if len(values) != 12:
        raise ValueError("expected length-12 list")
    return [int(values[k]) + int(values[k + 6]) for k in range(6)]


def f_2m_from_opposite_sums(s: list[int], m: int) -> complex:
    # n=12: F_{2m}(v) = Σ_{x=0}^{11} v(x) e^{i 4π m x / 12}
    # with opposite-pair reduction: F_{2m} = Σ_{k=0}^{5} (v(k)+v(k+6)) e^{i 2π m k / 6}.
    if len(s) != 6:
        raise ValueError("expected 6 opposite-pair sums")
    acc = 0.0 + 0.0j
    for k, sk in enumerate(s):
        ang = 2.0 * math.pi * float(m) * float(k) / 6.0
        acc += float(sk) * complex(math.cos(ang), math.sin(ang))
    return acc


def theta_base_from_real_fourier_defect(re_f: float) -> float:
    # For J(θ)=C + α*(Re(F_{2m})/12) cos(2θ) with α>0:
    # - if Re(F_{2m})>0: minima at cos(2θ)=-1 -> θ=π/2 (mod π)
    # - if Re(F_{2m})<0: minima at cos(2θ)=+1 -> θ=0 (mod π)
    if re_f == 0.0:
        raise ValueError("Re(F_{2m})=0: objective is constant on the O(2) orbit")
    return (math.pi / 2.0) if re_f > 0.0 else 0.0


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    alpha_geo_obj = load_json(IN_ALPHA_GEO)
    sigma_int_obj = load_json(IN_SIGMA_INT)

    sigma_int_val = sigma_int_obj.get("value")
    if not is_int_sign(sigma_int_val):
        raise SystemExit(
            json.dumps(
                {
                    "status": "INVALID_SIGMA_INT",
                    "expected": "sigma_int_strict_derived_v1.value ∈ {+1,-1}",
                    "actual": sigma_int_val,
                },
                ensure_ascii=True,
            )
        )

    alpha_geo_symbolic = str(((alpha_geo_obj.get("definition") or {}).get("computation")) or "4 ln(2)")
    r_ord_state = ensure_r_ord_present(alpha_geo_symbolic="4 ln(2)")
    r_ord = load_json(IN_R_ORD)

    ord_by_x = r_ord.get("ord_z12_by_x")
    if not (isinstance(ord_by_x, list) and len(ord_by_x) == 12 and all(isinstance(v, int) for v in ord_by_x)):
        raise SystemExit(
            json.dumps(
                {"status": "INVALID_R_ORD", "expected": "ord_z12_by_x: length-12 int list", "actual": ord_by_x},
                ensure_ascii=True,
            )
        )

    s = opposite_pair_sums([int(v) for v in ord_by_x])
    f2 = f_2m_from_opposite_sums(s, m=1)
    f4 = f_2m_from_opposite_sums(s, m=2)

    # These are provably real in N480/N488; we enforce numerical sanity only.
    f2_re, f2_im = float(f2.real), float(f2.imag)
    f4_re, f4_im = float(f4.real), float(f4.imag)
    if abs(f2_im) > 1e-9 or abs(f4_im) > 1e-9:
        raise SystemExit(
            json.dumps(
                {"status": "UNEXPECTED_NONREAL_FOURIER_DEFECT", "F2": {"re": f2_re, "im": f2_im}, "F4": {"re": f4_re, "im": f4_im}},
                ensure_ascii=True,
            )
        )

    theta_1_base = theta_base_from_real_fourier_defect(f2_re)  # expected π/2
    theta_2_base = theta_base_from_real_fourier_defect(f4_re)  # expected 0

    # Integrate the existing sigma-int Z2 sign convention from the export-map object (F311):
    # sigma_int=-1 corresponds to u_1 -> -u_1, i.e. theta_1 -> theta_1 + π.
    theta_1 = float(theta_1_base + (math.pi if int(sigma_int_val) == -1 else 0.0))
    theta_2 = float(theta_2_base)

    q2190 = load_json(IN_QW2190)
    n = int((q2190.get("mode_mapping") or {}).get("n_octaves") or 12)
    if n != 12:
        raise SystemExit(json.dumps({"status": "UNEXPECTED_N", "expected": 12, "actual": n}, ensure_ascii=True))

    fourier = real_fourier_basis(n)
    c1, s1 = fourier["c1"], fourier["s1"]
    c2, s2 = fourier["c2"], fourier["s2"]

    u1 = math.cos(theta_1) * c1 + math.sin(theta_1) * s1
    u2 = math.cos(theta_2) * c2 + math.sin(theta_2) * s2

    b = np.column_stack([u1, u2])

    theta_pair = {
        "object": "ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1",
        "status": "actual_exported_strict_core_theta_pair_source__slot_free_construction_class__no_false_pass",
        "as_of": "2026-03-15",
        "intent": (
            "Export one slot-free strict-core theta-pair source for the sigma-int corridor: "
            "sigma_int_strict_derived_v1 -> (theta_1, theta_2), derived from the strict Shannon element-order reference "
            "cross-entropy objective class (F446/N480) extended to pair2 (N488), with no eps/delta_d selector slots and no theta inputs."
        ),
        "typed_map_shape": "ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1 : sigma_int_strict_derived_v1 -> (theta_1, theta_2)",
        "construction_class": {
            "id": "T162_slot_free_sigma_int_to_theta_v1",
            "reference_distribution": str(IN_R_ORD.relative_to(REPO)),
            "objective_family": "J_{m}(θ) = -Σ_x p_{m,θ}(x) log r_ord(x)  (pair_m squared-amplitude distribution)",
            "theorem_refs": ["N479", "N480", "N488"],
            "slot_free": True,
            "no_theta_inputs": True,
            "no_eps_delta_d": True,
        },
        "inputs": {
            "sigma_int": {"object": "sigma_int_strict_derived_v1", "value": int(sigma_int_val)},
            "alpha_geo": {"object": "alpha_geo_strict_derived_v1", "value_symbolic": str(alpha_geo_obj.get("value"))},
        },
        "derived_fourier_defects_on_ord": {
            "ord_z12_by_x": [int(v) for v in ord_by_x],
            "opposite_pair_sums_Sk": [int(v) for v in s],
            "F2_ord_pair1": {"re": f2_re, "im": f2_im},
            "F4_ord_pair2": {"re": f4_re, "im": f4_im},
        },
        "outputs": {
            "pair1": {
                "m": 1,
                "theta_1_base_mod_pi": theta_1_base,
                "theta_1": theta_1,
                "sigma_int_sign_convention": "theta_1 := theta_1_base + π if sigma_int=-1 (matches F311)",
                "u_1": [float(x) for x in u1.tolist()],
            },
            "pair2": {
                "m": 2,
                "theta_2_base_mod_pi": theta_2_base,
                "theta_2": theta_2,
                "u_2": [float(x) for x in u2.tolist()],
            },
        },
        "audits": {
            "u1_u2_orthonormal_residual_vs_I2": orthonormal_residual(b),
            "u1_norm": float(np.linalg.norm(u1)),
            "u2_norm": float(np.linalg.norm(u2)),
            "u1_dot_u2": float(u1 @ u2),
        },
        "hard_limits": [
            "Does not discharge strict-derived eps selection (T160) nor strict-derived delta_d selection (T161).",
            "Does not claim admissible S_sel_int nor strict-core selector closure.",
            "Does not claim global QW-2191 discharge beyond the declared pair1/pair2 scope.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F451",
        "status": "F451_EXECUTED_CURRENT_STRICT_SIGMA_INT_SLOT_FREE_THETA_PAIR_SOURCE_PACKET_NO_FALSE_PASS",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "sigma_int_value": int(sigma_int_val),
            "alpha_geo_symbolic": str(alpha_geo_obj.get("value")),
            "r_ord_state": r_ord_state,
        },
        "outputs": {"theta_1": theta_1, "theta_2": theta_2, "u1_u2_orthonormal_residual": theta_pair["audits"]["u1_u2_orthonormal_residual_vs_I2"]},
        "no_false_pass": True,
    }

    OUT_THETA_PAIR.write_text(json.dumps(theta_pair, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_THETA_PAIR_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_THETA_PAIR_SUMMARY)


if __name__ == "__main__":
    main()

