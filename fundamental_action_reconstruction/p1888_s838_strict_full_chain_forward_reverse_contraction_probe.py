#!/usr/bin/env python3
"""P1888 S838 strict full-chain forward/reverse contraction probe."""
from __future__ import annotations
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1887 = load("p1887_s837_strict_nontrivial_transport_branch_probe.json")
    p1880 = load("p1880_s830_strict_full_lagrangian_term_registry_to_eom_transport_probe.json")

    beta, eta, omega, phi = sp.symbols("beta eta omega phi", positive=True, real=True)
    lam, y, m2 = sp.symbols("lam y m2", positive=True, real=True)

    c0 = sp.cos(omega + phi) / (beta + 1)
    c1 = (-beta * eta * sp.cos(omega + phi) - omega * (beta + 1) * sp.sin(omega + phi)) / (beta + 1) ** 2

    m2_eff = m2 * (1 + c0)
    lam_eff = lam * (1 + c1**2)
    y_eff = y * (1 + c0 / 2)

    # Reverse contraction placeholders: admissible corridor widths
    dm2, dlam, dy = sp.symbols("dm2 dlam dy", nonnegative=True, real=True)
    admissible_corridor = {
        "m2_eff_window": sp.Interval(m2_eff - dm2, m2_eff + dm2),
        "lam_eff_window": sp.Interval(lam_eff - dlam, lam_eff + dlam),
        "y_eff_window": sp.Interval(y_eff - dy, y_eff + dy),
    }

    qw2049 = {beta: 1.0, eta: 1.8, omega: 0.18575, phi: 0.16250}

    out = {
        "packet_id": "P1888",
        "stage_id": "S838",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1887_present": "nontrivial_transport_branch" in p1887,
            "p1880_present": "full_lagrangian_term_registry" in p1880,
        },
        "full_chain_forward": [
            "K_strict",
            "effective_coefficients",
            "full_non_skeleton_L_SM_plus_L_GR",
            "covariant_EOM",
            "one_loop_R/U checks",
            "FRW<->BianchiI_nontrivial_transport",
        ],
        "full_chain_reverse": [
            "R/U/BI_witness_solution",
            "term_registry_consistency",
            "EOM_residual_stability",
            "admissible_kernel_corridor",
        ],
        "effective_coefficients": {
            "m2_eff": str(m2_eff),
            "lam_eff": str(lam_eff),
            "y_eff": str(y_eff),
        },
        "reverse_contraction_contract": {
            "corridor_parameters": ["dm2", "dlam", "dy"],
            "admissible_windows": {k: str(v) for k, v in admissible_corridor.items()},
            "rule": "only witness-compatible coefficient windows are admissible for strict kernel continuation",
        },
        "qw2049_trace": {
            "m2_eff_over_m2": str(sp.N((m2_eff / m2).subs(qw2049), 12)),
            "lam_eff_over_lam": str(sp.N((lam_eff / lam).subs(qw2049), 12)),
            "y_eff_over_y": str(sp.N((y_eff / y).subs(qw2049), 12)),
        },
        "strict_core_closure_missing_items": {
            "renormalization": "Need solved one-loop finite/pole data to set corridor widths.",
            "unitarity": "Need solved Cutkosky defects on nontrivial transport branch.",
            "background_independence": "Need theorem-grade FRW<->BianchiI witness to legitimize reverse contraction.",
            "selector_obstruction": "QW-2191 remains open.",
        },
        "false_pass_guard": "Reverse contraction contract is not a solved strict-core closure theorem.",
        "next_honest_step": "Fill corridor widths from computed joint witness solution and test stability of admissible windows under branch perturbations.",
        "lay_explanation": "To krok spinający całość: z wyników testów kwantowych chcemy wrócić do tego, jakie ustawienia kernela są naprawdę dozwolone.",
    }

    path = GEN / "p1888_s838_strict_full_chain_forward_reverse_contraction_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
