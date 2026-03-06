#!/usr/bin/env python3
import json
from pathlib import Path


def main() -> int:
    payload = {
        "step": "H30",
        "status": "PASS_PARTIAL_KERNEL_INVARIANT_CANDIDATE_ONLY",
        "frontier": "H30_B1 := orientation_psi0 is a deterministic kernel-invariant anchor candidate derived from (phi,omega), but it is not yet exported or justified as the strict-core selector datum for theta_i and does not by itself discharge the residual O(2) degeneracy",
        "inputs": {
            "psi0_formula": "mod(0.5*phi + 0.8*omega, 2*pi)",
            "sources": [
                "QW_1952_INFORMATION_CHANNEL_DEDEGENERACY_OPERATOR.py",
                "QW_1953_TWO_STATE_INTERNAL_OBSERVER.py",
                "QW_1956_TWO_STATE_OBSERVER_WITH_REPAIRED_OPERATOR.py",
            ],
        },
        "classification": {
            "deterministic_from_kernel_invariants": True,
            "free_heuristic_parameter": False,
            "strict_core_selector_export": False,
            "residual_O2_discharge": False,
        },
        "hard_limits": {
            "theorem_level_pass": False,
            "full_closure_pass": False,
            "psi0_equals_selector_claim": False,
            "theta_export_claim": False,
            "qw2191_discharged": False,
        },
    }

    out_dir = Path(__file__).resolve().parent / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)
    for name in [
        "h30_kernel_invariant_psi0_anchor_audit.json",
        "h30_kernel_invariant_psi0_anchor_audit_summary.json",
    ]:
        (out_dir / name).write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
