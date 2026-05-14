from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

TOL = 1e-8


def coeffs_from_params(omega: float, phi: float, beta: float, eta: float) -> tuple[float, float, float]:
    lam = abs(math.cos(phi)) / (1.0 + 0.3 * beta)
    kap = (omega * omega + 0.5 * eta) / (0.2 + abs(math.cos(phi)))
    eps = (beta * eta + 0.1 * omega) / (1.0 + abs(math.sin(phi)))
    return (lam, kap, eps)


def xi_phase_curv(omega: float, phi: float, beta: float, eta: float) -> float:
    return omega * eta + math.sin(phi)


def l2(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def solve_3x3(a: list[list[float]], b: list[float]) -> list[float]:
    m = [row[:] + [rhs] for row, rhs in zip(a, b)]
    n = 3
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(m[r][i]))
        m[i], m[piv] = m[piv], m[i]
        p = m[i][i]
        if abs(p) < 1e-14:
            raise ValueError("Singular Jacobian in Newton solve")
        for j in range(i, n + 1):
            m[i][j] /= p
        for r in range(n):
            if r == i:
                continue
            fac = m[r][i]
            for j in range(i, n + 1):
                m[r][j] -= fac * m[i][j]
    return [m[i][n] for i in range(n)]


def fit_point_with_fixed_eta(target: tuple[float, float, float], eta: float) -> tuple[float, float, float]:
    x = [0.30, 0.15, 1.10]  # omega, phi, beta
    h = 1e-6
    for _ in range(40):
        f = coeffs_from_params(x[0], x[1], x[2], eta)
        r = [f[i] - target[i] for i in range(3)]
        if max(abs(v) for v in r) < 1e-12:
            break
        cols: list[list[float]] = []
        for j in range(3):
            xp = x[:]
            xp[j] += h
            fp = coeffs_from_params(xp[0], xp[1], xp[2], eta)
            cols.append([(fp[i] - f[i]) / h for i in range(3)])
        jac = [[cols[c][r_idx] for c in range(3)] for r_idx in range(3)]
        dx = solve_3x3(jac, [-v for v in r])
        x = [x[i] + dx[i] for i in range(3)]
    return (x[0], x[1], x[2])


def main() -> None:
    root = Path(__file__).resolve().parent
    generated = root / "generated"
    generated.mkdir(parents=True, exist_ok=True)

    p0 = (0.18575, 0.16250, 1.0, 1.8)
    c0 = coeffs_from_params(*p0)
    xi0 = xi_phase_curv(*p0)

    eta1 = 1.65
    omega1, phi1, beta1 = fit_point_with_fixed_eta(c0, eta1)
    p1 = (omega1, phi1, beta1, eta1)
    c1 = coeffs_from_params(*p1)
    xi1 = xi_phase_curv(*p1)

    coeff_distance = l2(c0, c1)
    params_distance = math.sqrt(sum((a - b) ** 2 for a, b in zip(p0, p1)))

    w1567a_pass = coeff_distance < TOL and params_distance > 1e-2
    xi_gap = abs(xi0 - xi1)
    w1567b_pass = xi_gap > 1e-2

    status = (
        "PASS_W1567A_AND_PASS_W1567B_CANDIDATE"
        if (w1567a_pass and w1567b_pass)
        else "FAIL_W1567A_ORBIT_NOT_WITNESSED"
    )

    summary = {
        "checkpoint": "P1568_S518",
        "date_utc": "2026-05-14",
        "status": status,
        "strict_only": True,
        "legacy_bridge_used": False,
        "qw2191_closed": False,
        "toe_closed": False,
        "base_params": {"omega": p0[0], "phi": p0[1], "beta": p0[2], "eta": p0[3]},
        "alt_params": {"omega": p1[0], "phi": p1[1], "beta": p1[2], "eta": p1[3]},
        "base_coeffs": {"lambda_sm_eff": c0[0], "kappa_gr_eff": c0[1], "epsilon_mix_eff": c0[2]},
        "alt_coeffs": {"lambda_sm_eff": c1[0], "kappa_gr_eff": c1[1], "epsilon_mix_eff": c1[2]},
        "coeff_bundle_l2_distance": coeff_distance,
        "param_space_l2_distance": params_distance,
        "W1567A_pass": w1567a_pass,
        "xi_phase_curv_base": xi0,
        "xi_phase_curv_alt": xi1,
        "xi_phase_curv_gap": xi_gap,
        "W1567B_pass": w1567b_pass,
        "tolerances": {"coeff_bundle_l2": TOL, "xi_gap_min": 1e-2},
        "missing_exports_witnesses_theorems": [
            "E1568C_analytic_gradient_hessian_export_for_xi",
            "T1568D_local_uniqueness_theorem_for_4observable_map",
            "T1568E_inverse_stability_envelope_theorem",
        ],
        "next_honest_step": "P1569_formalize_local_uniqueness_theorem_T1568D",
        "lay_summary": "Znaleziono dwie różne konfiguracje strict o tym samym bundle 3x1; 4-ty sygnał xi_phase_curv je rozdziela.",
    }

    summary["audit_digest"] = hashlib.sha256(json.dumps(summary, sort_keys=True).encode("utf-8")).hexdigest()

    out = generated / "p1568_s518_strict_underdetermined_orbit_witness_and_4th_observable_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1568] wrote {out} status={status}")


if __name__ == "__main__":
    main()
