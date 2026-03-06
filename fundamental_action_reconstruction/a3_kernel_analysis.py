#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

AS_OF = "2026-03-06"

summary = {
    "program_id": "fundamental_action_reconstruction",
    "as_of": AS_OF,
    "status": "A3_MINIMAL_BRANCH_OPERATOR_SPLIT_COMPLETED_NO_FULL_CLOSURE_CLAIM",
    "anti_overclaim": {
        "global_stability_claim": False,
        "full_bosonic_positivity_claim": False,
        "fermionic_kernel_closed_claim": False,
        "lorentzian_unitarity_claim": False,
    },
    "input_branch": "single_foundation_minimal_matching_branch",
    "fluctuations": {
        "physical_doublet": ["delta_rho", "delta_phi"],
        "orthogonal_shape_sector": "delta n_perp^A after projection away from moduli",
        "gauge_sector": "inactive on the executed branch; only deferred as a future projection layer",
    },
    "operator": {
        "physical_operator": "O_phys = - d/dr [ K_2(r) d/dr ] + M_2(r)",
        "kinetic_block": [
            ["r^(d-1) K(rho_bg,phi_bg)", "0"],
            ["0", "r^(d-1)"],
        ],
        "mass_block": [
            "r^(d-1) Hess_(rho,phi) V_eff evaluated on the background",
            "- r^(d-3) Hess_(rho,phi) T_top evaluated on the background",
            "+ connection/background-gradient corrections inherited from K(rho,phi)",
        ],
    },
    "mode_split": {
        "zero_modes": [
            "translational moduli if the background is localized",
            "internal orientation moduli if n^A has a continuous moduli manifold",
            "scale mode only on a critical/BPS-like branch; not assumed generically",
        ],
        "gauge_modes": [
            "none in the currently executed gauge-off variable set",
            "must be projected out once A_mu^I is activated on later branches",
        ],
        "physical_modes": [
            "delta_rho amplitude channel",
            "delta_phi order-parameter channel",
            "mixed amplitude/order channel from M_2(r)",
            "orthogonal shape modes after zero-mode projection",
        ],
    },
    "stability_criteria": {
        "bosonic_euclidean": [
            "project out zero modes first",
            "project out gauge modes once the gauge sector is active",
            "check coercivity/positivity only on the remaining physical subspace",
            "exclude tachyonic instability in the declared scope",
        ],
        "fermionic": [
            "do not require positivity of the Dirac operator itself",
            "future checks must use self-adjointness / ellipticity / positivity of D^dagger D / reflection positivity",
        ],
        "lorentzian": [
            "future checks must use well-posedness",
            "hyperbolicity",
            "bounded-below Hamiltonian",
            "absence of ghosts",
        ],
    },
    "open_obligations": [
        "compute coercivity of O_phys for a concrete background profile",
        "decide whether orthogonal shape modes are lifted or remain moduli",
        "extend the split once the gauge sector is activated",
        "carry the physical operator into a coarse-graining step in A4",
    ],
    "next_steps": ["A4", "A5"],
}

root = Path(__file__).resolve().parent
out = root / "generated" / "a3_kernel_analysis_summary.json"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
print(f"wrote {out}")
