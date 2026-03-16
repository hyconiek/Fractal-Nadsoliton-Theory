from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

root = Path(__file__).resolve().parent
repo = root.parent
generated = root / "generated"

R12_SUMMARY = generated / "r12_explicit_canonical_psi_block_export_packet_with_kernel_mixing_for_host_route_summary.json"
M_CONTROL_OBJECT = generated / "m_control_canonical_psi_declared_coefficient_filled_v1.json"
P476_SUMMARY = (
    generated
    / "p476_current_strict_r11_r12_declared_control_pullback_m_control_coefficient_filled_export_probe_summary.json"
)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _p(path: Path) -> str:
    try:
        return str(path.relative_to(repo))
    except ValueError:
        return str(path)


def main() -> None:
    generated.mkdir(exist_ok=True)

    coefficient_filled_h_present = False
    coefficient_filled_h_source: str | None = None
    if R12_SUMMARY.exists():
        r12s = _read_json(R12_SUMMARY)
        coefficient_filled_h_present = str(r12s.get("status", "")).startswith(
            "PASS_PARTIAL_EXPLICIT_CANONICAL_PSI_BLOCK_EXPORT"
        )
        if coefficient_filled_h_present:
            coefficient_filled_h_source = _p(R12_SUMMARY)

    coefficient_filled_m_present = bool(M_CONTROL_OBJECT.exists())
    coefficient_filled_m_source: str | None = _p(M_CONTROL_OBJECT) if coefficient_filled_m_present else None
    if coefficient_filled_m_present and P476_SUMMARY.exists():
        # Keep a stable pointer to the producing stage when available.
        coefficient_filled_m_source = _p(P476_SUMMARY)

    residual_blockers: dict[str, str] = {
        "C15_B2": "no_explicit_restriction_from_M_control_to_the_candidate_orientation_slice",
        "C9_B2": "no_explicit_restriction_from_that_Psi_sector_quadratic_carrier_to_the_candidate_orientation_slice",
    }
    if not coefficient_filled_h_present:
        residual_blockers = {
            "C15_B1": "no_explicit_coefficient_filled_canonical_Psi_x_Psi_block_H_PsiPsi_for_evaluating_the_control_pullback",
            **residual_blockers,
        }

    payload = {
        "stage": "C15",
        "status": "C15_EXECUTED_CONTROL_ONLY_PULLBACK_PACKET_NO_FALSE_PASS",
        "as_of": datetime.now(timezone.utc).date().isoformat(),
        "goal": (
            "Track the strict control-only pullback assembly target "
            "M_control = T_control^T H_PsiPsi T_control and whether coefficient-filled artifacts exist "
            "for H_PsiPsi and M_control in declared scope."
        ),
        "inputs": {
            "strict_admissible": [
                "QW-2190",
                "QW-2164",
                "QW-2166",
                "QW-2180",
                "C12",
                "C14",
                "A10",
                "R11",
                "R12",
            ]
        },
        "formal_objects": {
            "T_control_shape": [12, 4],
            "H_PsiPsi_shape": [12, 12],
            "M_control_shape": [4, 4],
            "control_basis": ["c1", "s1", "c2", "s2"],
            "assembly_formula": "M_control = T_control^T H_PsiPsi T_control",
        },
        "result": {
            "control_only_pullback_formula_present": "yes",
            "coefficient_filled_H_PsiPsi_present": "yes" if coefficient_filled_h_present else "not_shown",
            "coefficient_filled_H_PsiPsi_source": coefficient_filled_h_source,
            "coefficient_filled_M_control_present": "yes" if coefficient_filled_m_present else "not_shown",
            "coefficient_filled_M_control_source": coefficient_filled_m_source,
            "orientation_slice_restriction_present": "not_shown",
        },
        "residual_blockers": residual_blockers,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_physical_canonicalization_claim",
            "no_host_matching_claim",
            "no_orientation_slice_restriction_claim",
            "no_qw2191_discharge_claim",
        ],
        "next_step": "C26",
        "no_false_pass": True,
    }

    out = generated / "c15_control_only_pullback_submatrix_packet_summary.json"
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(out)


if __name__ == "__main__":
    main()
