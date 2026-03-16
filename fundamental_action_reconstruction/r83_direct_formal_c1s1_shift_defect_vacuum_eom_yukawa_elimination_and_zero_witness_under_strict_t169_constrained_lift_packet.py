#!/usr/bin/env python3
from __future__ import annotations

import json
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

# Inputs
R14 = GENERATED / "r14_explicit_frozen_kernel_specialization_packet_for_host_matching_route.json"
R21 = GENERATED / "r21_explicit_pair1_c1s1_shift_defect_polynomial_packet_for_host_route.json"
F447 = GENERATED / "f447_current_strict_t169_qw2122_scalar_to_t168_per_site_value_provider_strict_derived_v1.json"

N474 = ROOT / "N474_CURRENT_FIRST_STRICT_VACUUM_EOM_ELIMINATES_YUKAWA_FROM_CANONICAL_LOCAL_DIAGONAL_RESIDUAL_THEOREM.md"
N483 = ROOT / "N483_CURRENT_FIRST_STRICT_T169_POWERLAW_ELEMENT_ORDER_CONSTRAINED_LIFT_EXISTENCE_UNIQUENESS_THEOREM.md"

# Outputs
OUT_JSON = (
    GENERATED
    / "r83_direct_formal_c1s1_shift_defect_vacuum_eom_yukawa_elimination_and_zero_witness_under_strict_t169_constrained_lift_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "r83_direct_formal_c1s1_shift_defect_vacuum_eom_yukawa_elimination_and_zero_witness_under_strict_t169_constrained_lift_packet_summary.json"
)


def load_json_decimal(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"), parse_float=Decimal)


def inv_z12(i: int) -> int:
    return (-i) % 12


def parse_slot_index(slot: str) -> int:
    if not slot.startswith("psi"):
        raise ValueError(f"unexpected carrier slot label: {slot}")
    return int(slot[3:])


def max_abs(xs: list[Decimal]) -> Decimal:
    if not xs:
        return Decimal(0)
    return max(abs(x) for x in xs)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)
    getcontext().prec = 80

    missing: list[str] = []
    for dep in (R14, R21, F447):
        if not dep.exists():
            missing.append(str(dep.relative_to(REPO)))
    if missing:
        payload = {
            "stage": "R83",
            "status": "FAIL_MISSING_DEPENDENCY_FILES",
            "missing_dependency_files": missing,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {"stage": "R83", "status": payload["status"], "shift_defect": None, "no_false_pass": True},
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_SUMMARY)
        return

    r14 = load_json_decimal(R14)
    r21 = load_json_decimal(R21)
    f447 = load_json_decimal(F447)

    host_kernel_rows = r14.get("host_kernel_rows") or []
    if not (isinstance(host_kernel_rows, list) and len(host_kernel_rows) == 12 and all(len(r) == 12 for r in host_kernel_rows)):
        raise SystemExit(
            json.dumps(
                {
                    "stage": "R83",
                    "status": "FAIL_INVALID_R14_HOST_KERNEL_ROWS_SHAPE",
                    "r14_path": str(R14.relative_to(REPO)),
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )
    K_total: list[list[Decimal]] = host_kernel_rows

    pos_slots = (r21.get("pair1_c1s1_shift_defect_packet") or {}).get("positive_support_slots") or []
    neg_slots = (r21.get("pair1_c1s1_shift_defect_packet") or {}).get("negative_support_slots") or []

    pos = [parse_slot_index(s) for s in pos_slots]
    neg = [parse_slot_index(s) for s in neg_slots]

    if not (len(pos) == 4 and len(neg) == 4):
        raise SystemExit(
            json.dumps(
                {
                    "stage": "R83",
                    "status": "FAIL_UNEXPECTED_C1S1_SUPPORT_CARDINALITY",
                    "positive_support_slots": pos_slots,
                    "negative_support_slots": neg_slots,
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    pos_inv = [inv_z12(i) for i in pos]
    support_inv_matches = sorted(pos_inv) == sorted(neg)

    vpsi = [Decimal(x) for x in (f447.get("vpsi") or [])]
    g4 = [Decimal(x) for x in (f447.get("g4") or [])]
    g6 = [Decimal(x) for x in (f447.get("g6") or [])]
    values = f447.get("values") or {}
    m0_sq = Decimal(values.get("m0_sq"))

    if not (len(vpsi) == 12 and len(g4) == 12 and len(g6) == 12):
        raise SystemExit(
            json.dumps(
                {
                    "stage": "R83",
                    "status": "FAIL_INVALID_F447_ARRAY_SHAPES",
                    "vpsi_len": len(vpsi),
                    "g4_len": len(g4),
                    "g6_len": len(g6),
                    "f447_path": str(F447.relative_to(REPO)),
                    "no_false_pass": True,
                },
                ensure_ascii=True,
            )
        )

    vpsi_nonzero_all = all(x != 0 for x in vpsi)
    min_abs_vpsi = min(abs(x) for x in vpsi) if vpsi else Decimal(0)

    max_abs_vpsi_inv_diff = max_abs([vpsi[i] - vpsi[inv_z12(i)] for i in range(12)])
    max_abs_ktotal_inv_diff = max_abs(
        [K_total[i][j] - K_total[inv_z12(i)][inv_z12(j)] for i in range(12) for j in range(12)]
    )

    g4_uniform = len(set(g4)) == 1
    g6_zero = all(x == 0 for x in g6)

    def mixing_ratio(k: int) -> Decimal:
        s = Decimal(0)
        for j in range(12):
            if j == k:
                continue
            s += K_total[k][j] * (vpsi[j] / vpsi[k])
        return s

    def d_local_residual_yukawa_elim(k: int) -> Decimal:
        return (
            -mixing_ratio(k)
            + Decimal(2) * g4[k] * (vpsi[k] ** 2)
            + Decimal(4) * g6[k] * (vpsi[k] ** 4)
            - m0_sq
        )

    d = [d_local_residual_yukawa_elim(k) for k in range(12)]
    max_abs_d_inv_diff = max_abs([d[i] - d[inv_z12(i)] for i in range(12)])

    shift_defect = sum(d[i] for i in pos) - sum(d[i] for i in neg)

    checks = [
        {
            "id": "n474_theorem_exists",
            "actual": bool(N474.exists()),
            "expected": True,
            "meaning": "N474 exports the conditional vacuum-EoM Yukawa elimination tool for the diagonal residual entry d_k",
        },
        {
            "id": "n483_theorem_exists",
            "actual": bool(N483.exists()),
            "expected": True,
            "meaning": "N483 exports the strict T169 constrained lift rule used by F447",
        },
        {
            "id": "support_sets_match_under_inversion_i_to_minus_i_mod_12",
            "actual": support_inv_matches,
            "expected": True,
            "meaning": "the c1s1 support positive and negative index-sets are related by inversion i ↦ -i (mod 12)",
        },
        {
            "id": "f447_vpsi_nonzero_all_sites",
            "actual": bool(vpsi_nonzero_all),
            "expected": True,
            "meaning": "N474 applies sitewise only when vpsi_k ≠ 0; the exported T169 instance should satisfy this on all sites used here",
        },
        {
            "id": "f447_vpsi_inversion_invariant_exported_instance",
            "actual": bool(max_abs_vpsi_inv_diff == 0),
            "expected": True,
            "actual_value": str(max_abs_vpsi_inv_diff),
            "meaning": "the exported T169 constrained-lift vacuum vector vpsi is inversion-invariant: vpsi[i] = vpsi[-i mod 12]",
        },
        {
            "id": "r14_ktotal_inversion_invariant",
            "actual": bool(max_abs_ktotal_inv_diff == 0),
            "expected": True,
            "actual_value": str(max_abs_ktotal_inv_diff),
            "meaning": "the frozen K_total carrier is inversion-invariant: K_total[i,j] = K_total[-i mod 12, -j mod 12]",
        },
        {
            "id": "f447_g4_uniform",
            "actual": bool(g4_uniform),
            "expected": True,
            "meaning": "under T169 the quartic family g4_psi_i is exported as uniform (F447)",
        },
        {
            "id": "f447_g6_uniform_zero",
            "actual": bool(g6_zero),
            "expected": True,
            "meaning": "under T169 the sextic family g6_psi_i is exported as identically zero (F447)",
        },
        {
            "id": "yukawa_eliminated_shift_defect_zero_on_exported_instance",
            "actual": bool(shift_defect == 0),
            "expected": True,
            "actual_value": str(shift_defect),
            "meaning": (
                "under the N474 elimination form, the evaluated c1s1 shift defect Σ_pos d_k - Σ_neg d_k is zero on the exported instance "
                "(R14 + F447), hence the c1s1 defect is closed on this declared instance under the explicit vacuum-EoM premise scope"
            ),
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    status = "PASS_R83_SHIFT_DEFECT_ZERO_WITNESS_UNDER_N474_YUKAWA_ELIMINATION_AND_T169_CONSTRAINED_LIFT"
    if not all(bool(c["pass"]) for c in checks):
        status = "FAIL_R83_CHECKS_NOT_SATISFIED"

    artifact = {
        "stage": "R83",
        "goal": (
            "close_pair1_c1s1_shift_defect_under_vacuum_eom_yukawa_elimination_by_evaluating_the_yukawa_free_diagonal_residual_form "
            "on_the_exported_strict_T169_constrained_lift_instance_and_frozen_K_total_specialization"
        ),
        "dependencies": {
            "R14": str(R14.relative_to(REPO)),
            "R21": str(R21.relative_to(REPO)),
            "F447": str(F447.relative_to(REPO)),
            "N474": str(N474.relative_to(REPO)),
            "N483": str(N483.relative_to(REPO)),
        },
        "premise_scope_note": (
            "R83 uses N474, which is conditional on constant-vacuum stationarity and requires vpsi_k ≠ 0. "
            "This packet is an exported-instance witness only and does not claim a global theorem."
        ),
        "c1s1_support": {
            "positive_support_slots": pos_slots,
            "negative_support_slots": neg_slots,
            "positive_indices": pos,
            "negative_indices": neg,
            "positive_indices_inverted": pos_inv,
            "support_inversion_matches": support_inv_matches,
        },
        "instance_invariance_audits": {
            "min_abs_vpsi": str(min_abs_vpsi),
            "max_abs_vpsi_inversion_diff": str(max_abs_vpsi_inv_diff),
            "max_abs_ktotal_inversion_diff": str(max_abs_ktotal_inv_diff),
            "max_abs_d_local_residual_inversion_diff_under_elimination_form": str(max_abs_d_inv_diff),
        },
        "yukawa_eliminated_form_used": {
            "d_k_definition": (
                "d_k := -Σ_{j≠k} K_total[k,j] * (vpsi_j/vpsi_k) + 2*g4_k*vpsi_k^2 + 4*g6_k*vpsi_k^4 - m0^2  (from N474 + R15)"
            ),
            "note": "This form is Yukawa-free: gY_k and vphi do not appear.",
        },
        "result": {
            "status": status,
            "shift_defect_value": str(shift_defect),
            "shift_defect_zero": bool(shift_defect == 0),
            "pair1_c1s1_shift_defect_zero_witness_present_under_elimination_form": bool(shift_defect == 0),
            "strict_core_promotion": False,
        },
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_strict_core_promotion",
            "does_not_claim_pair1_c1c1_or_s1s1_equations",
            "does_not_claim_QW2191_discharge",
            "does_not_claim_selector_closure",
            "does_not_claim_ToE_closure",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["result"]["status"],
        "shift_defect_value": artifact["result"]["shift_defect_value"],
        "shift_defect_zero": artifact["result"]["shift_defect_zero"],
        "pair1_c1s1_shift_defect_zero_witness_present_under_elimination_form": artifact["result"][
            "pair1_c1s1_shift_defect_zero_witness_present_under_elimination_form"
        ],
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
