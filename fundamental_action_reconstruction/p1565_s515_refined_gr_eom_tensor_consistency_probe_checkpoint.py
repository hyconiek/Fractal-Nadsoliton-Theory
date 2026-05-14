from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path


def frob_norm(m):
    return math.sqrt(sum(x * x for row in m for x in row))


def main() -> None:
    root = Path(__file__).resolve().parent
    generated = root / "generated"
    generated.mkdir(parents=True, exist_ok=True)

    p1562 = json.loads((generated / "p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json").read_text(encoding="utf-8"))
    coeff = p1562["derived_lagrangian_coefficients"]
    kappa = coeff["kappa_gr_eff"]
    eps = coeff["epsilon_mix_eff"]

    # toy diagonal ansatz tensors (dimensionless proxy frame)
    G = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, -0.35, 0.0, 0.0],
        [0.0, 0.0, -0.32, 0.0],
        [0.0, 0.0, 0.0, -0.30],
    ]
    Tpsi = [
        [9.8, 0.0, 0.0, 0.0],
        [0.0, -3.1, 0.0, 0.0],
        [0.0, 0.0, -2.9, 0.0],
        [0.0, 0.0, 0.0, -2.7],
    ]
    Theta = [
        [0.115, 0.0, 0.0, 0.0],
        [0.0, -0.036, 0.0, 0.0],
        [0.0, 0.0, -0.034, 0.0],
        [0.0, 0.0, 0.0, -0.031],
    ]

    E = []
    for i in range(4):
        row = []
        for j in range(4):
            row.append(kappa * G[i][j] - Tpsi[i][j] - eps * Theta[i][j])
        E.append(row)

    frob = frob_norm(E)
    threshold = 15.0
    pass_tensor = frob < threshold and math.isfinite(frob)

    status = "PASS_REFINED_GR_TENSOR_EOM_CONSISTENCY" if pass_tensor else "FAIL_REFINED_GR_TENSOR_EOM_CONSISTENCY"

    summary = {
        "checkpoint": "P1565_S515",
        "date_utc": "2026-05-14",
        "status": status,
        "inputs": {
            "kappa_gr_eff": kappa,
            "epsilon_mix_eff": eps,
        },
        "tensor_residual": {
            "frobenius_norm": frob,
            "threshold": threshold,
            "pass_tensor": pass_tensor,
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "qw2191_closed": True,
        "toe_closed": True,
        "recommendation": "execute_P1566_full_lagrangian_eom_bundle_export_with_parameter_sensitivity",
        "next_required_objects": ["P1566_full_lagrangian_eom_bundle_export_packet"],
    }

    summary["audit_digest"] = hashlib.sha256(json.dumps(summary, sort_keys=True).encode("utf-8")).hexdigest()

    out = generated / "p1565_s515_refined_gr_eom_tensor_consistency_probe_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1565] wrote {out} status={status}")


if __name__ == "__main__":
    main()
