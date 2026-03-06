#!/usr/bin/env python3
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated"
OUT.mkdir(exist_ok=True)

payload = {
    "step": "V3",
    "status": "PASS_PARTIAL_COMPETING_EXTENSION_HYPOTHESIS_ONLY",
    "date": "2026-03-06",
    "carrier": {
        "pair": "pair1",
        "space": "V_1 = span{c_1,s_1}",
        "basis": ["c_1", "s_1"],
    },
    "admissible_classes": {
        "isotropic": {
            "formula": "K_visc_iso^(1) = nu_iso * I_2",
            "requires_anchor": False,
            "can_break_O2": False,
        },
        "anchor_imported_anisotropic": {
            "formula": "K_visc_aniso^(1)(psi0) = R(psi0) diag(nu_parallel, nu_perp) R(-psi0)",
            "requires_anchor": True,
            "anchor_source": "external_candidate_only",
            "can_break_O2": True,
        },
    },
    "selector_smuggling_forbidden": [
        "theta_star_fixed_by_definition",
        "preferred_basis_direction_without_exported_anchor",
        "anisotropy_renamed_as_viscosity_without_explanation",
    ],
    "frontier": "V3_B1 := a minimal pair-level informational-viscosity operator can be written on V_1, but only in two admissible classes: isotropic (which cannot break O(2)) or anchor-imported anisotropic (which depends on an external anchor and therefore does not by itself solve selector closure)",
    "hard_limits": {
        "theorem_level_pass": False,
        "full_closure_pass": False,
        "viscosity_already_in_strict_core_as_selector": False,
        "isotropic_viscosity_breaks_o2": False,
        "anisotropic_class_self_generates_anchor": False,
    },
}

json_path = OUT / "v3_minimal_pair_level_viscosity_operator_audit.json"
summary_path = OUT / "v3_minimal_pair_level_viscosity_operator_audit_summary.json"
json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
summary_path.write_text(json.dumps({
    "step": payload["step"],
    "status": payload["status"],
    "carrier": payload["carrier"],
    "admissible_classes": payload["admissible_classes"],
    "frontier": payload["frontier"],
    "hard_limits": payload["hard_limits"],
}, indent=2, ensure_ascii=False) + "\n")
