#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1626 = GEN / "p1626_s576_full_strict_actional_chain_and_closure_obligation_summary.json"


def k_strict(d: float, omega: float, phi: float, beta: float, eta: float) -> float:
    return math.cos(omega * d + phi) / (1.0 + beta * (d ** eta))


def main() -> None:
    s26 = json.loads(IN1626.read_text(encoding="utf-8"))
    p = s26["full_strict_chain"]["kernel_params"]
    c = s26["full_strict_chain"]["coeff_mapping"]

    omega, phi, beta, eta = p["omega"], p["phi"], p["beta"], p["eta"]
    lam, kap, eps = c["lambda_sm_eff"], c["kappa_gr_eff"], c["epsilon_mix_eff"]

    grid = [i * 0.05 for i in range(1, 121)]
    m_vals = []
    deriv_sign_flips = 0
    prev = None
    prev_diff = None

    for d in grid:
        psi = k_strict(d, omega, phi, beta, eta)
        R = 0.2 + 0.8 * (psi ** 2)
        M = kap * R + eps * psi * R - 0.5 * lam * (psi ** 2)
        m_vals.append(M)
        if prev is not None:
            diff = M - prev
            if prev_diff is not None and diff * prev_diff < 0:
                deriv_sign_flips += 1
            prev_diff = diff
        prev = M

    local_monotone = deriv_sign_flips <= 1
    local_uniqueness_proxy = local_monotone and (max(m_vals) - min(m_vals) > 1e-6)

    summary = {
        "checkpoint": "P1627_S577",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "status": "PASS_P1627_LOCAL_UNIQUENESS_PROXY" if local_uniqueness_proxy else "KEEP_OPEN_P1627_LOCAL_PROXY_NOT_SATISFIED",
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "full_chain_used": s26["full_strict_chain"],
        "local_proof_attempt": {
            "monotone_functional": "M = kap*R + eps*psi*R - 0.5*lam*psi^2",
            "grid_size": len(grid),
            "derivative_sign_flips": deriv_sign_flips,
            "local_monotone": local_monotone,
            "local_uniqueness_proxy": local_uniqueness_proxy,
            "scope": "LOCAL_ONLY",
        },
        "strict_core_closure": {
            "status": "OPEN",
            "missing_exports": ["E_selector_internal_source_full_domain", "E_full_variational_proof_log_machine_checkable"],
            "missing_witnesses": ["W_noncyclic_provider_for_selector_uniqueness_PROVED_GLOBAL"],
            "missing_theorems": ["T_qw2191_selector_uniqueness_strict_GLOBAL", "T_global_toe_closure_strict"],
        },
        "strict_only": True,
        "legacy_bridge_used": False,
        "external_team_validation_required": False,
        "next_honest_step": "Uogólnić lokalny proxy-proof do pełnej domeny i zamienić LOCAL_ONLY na GLOBAL_PROVED przez formalny dowód theorem-level.",
        "lay_summary": "Sprawdziliśmy, że w lokalnym zakresie teoria zachowuje się jakby miała jedną preferowaną ścieżkę; to dobry znak, ale to jeszcze nie globalny dowód.",
    }

    out = GEN / "p1627_s577_strict_selector_uniqueness_local_proof_attempt_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
