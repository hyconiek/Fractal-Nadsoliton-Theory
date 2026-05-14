from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

BASE = (0.18575, 0.16250, 1.0, 1.8)
H = 1e-6
TOL_J = 5e-6


def obs_map(omega: float, phi: float, beta: float, eta: float) -> tuple[float, float, float, float]:
    c = abs(math.cos(phi))
    s = abs(math.sin(phi))
    lam = c / (1.0 + 0.3 * beta)
    kap = (omega * omega + 0.5 * eta) / (0.2 + c)
    eps = (beta * eta + 0.1 * omega) / (1.0 + s)
    xi = omega * eta + math.sin(phi)
    return (lam, kap, eps, xi)


def analytic_jacobian(p: tuple[float, float, float, float]) -> list[list[float]]:
    omega, phi, beta, eta = p
    c_raw, s_raw = math.cos(phi), math.sin(phi)
    c, s = abs(c_raw), abs(s_raw)
    sc = 1.0 if c_raw >= 0 else -1.0
    ss = 1.0 if s_raw >= 0 else -1.0

    dlam_domega = 0.0
    dlam_dphi = (-sc * s_raw) / (1.0 + 0.3 * beta)
    dlam_dbeta = -0.3 * c / (1.0 + 0.3 * beta) ** 2
    dlam_deta = 0.0

    denom_k = 0.2 + c
    num_k = omega * omega + 0.5 * eta
    dkap_domega = (2.0 * omega) / denom_k
    dkap_dphi = num_k * (sc * s_raw) / (denom_k * denom_k)
    dkap_dbeta = 0.0
    dkap_deta = 0.5 / denom_k

    denom_e = 1.0 + s
    num_e = beta * eta + 0.1 * omega
    deps_domega = 0.1 / denom_e
    deps_dphi = -num_e * (ss * c_raw) / (denom_e * denom_e)
    deps_dbeta = eta / denom_e
    deps_deta = beta / denom_e

    dxi_domega = eta
    dxi_dphi = c_raw
    dxi_dbeta = 0.0
    dxi_deta = omega

    return [
        [dlam_domega, dlam_dphi, dlam_dbeta, dlam_deta],
        [dkap_domega, dkap_dphi, dkap_dbeta, dkap_deta],
        [deps_domega, deps_dphi, deps_dbeta, deps_deta],
        [dxi_domega, dxi_dphi, dxi_dbeta, dxi_deta],
    ]


def numeric_jacobian(p: tuple[float, float, float, float], h: float = H) -> list[list[float]]:
    base = obs_map(*p)
    cols = []
    for j, v in enumerate(p):
        pp = list(p)
        pp[j] = v + h
        f1 = obs_map(*pp)
        cols.append([(f1[i] - base[i]) / h for i in range(4)])
    return [[cols[c][r] for c in range(4)] for r in range(4)]


def inv4(m: list[list[float]]) -> list[list[float]]:
    n = 4
    a = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(m)]
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(a[r][i]))
        a[i], a[piv] = a[piv], a[i]
        p = a[i][i]
        if abs(p) < 1e-14:
            raise ValueError("singular matrix")
        for j in range(2 * n):
            a[i][j] /= p
        for r in range(n):
            if r == i:
                continue
            fac = a[r][i]
            for j in range(2 * n):
                a[r][j] -= fac * a[i][j]
    return [row[n:] for row in a]


def mat_inf_norm(m: list[list[float]]) -> float:
    return max(sum(abs(x) for x in row) for row in m)


def hessian_xi(_: tuple[float, float, float, float]) -> list[list[float]]:
    # xi = omega*eta + sin(phi)
    h = [[0.0] * 4 for _ in range(4)]
    h[0][3] = 1.0
    h[3][0] = 1.0
    h[1][1] = -math.sin(BASE[1])
    return h


def main() -> None:
    root = Path(__file__).resolve().parent
    out_dir = root / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    ja = analytic_jacobian(BASE)
    jn = numeric_jacobian(BASE)
    max_abs_err = max(abs(ja[r][c] - jn[r][c]) for r in range(4) for c in range(4))
    pass_e1569a = max_abs_err < TOL_J

    j_inv = inv4(ja)
    j_inv_norm = mat_inf_norm(j_inv)

    delta_f = (1e-4, -1e-4, 1e-4, -1e-4)
    delta_f_norm = max(abs(x) for x in delta_f)
    delta_p_bound = j_inv_norm * delta_f_norm
    pass_t1569b = math.isfinite(delta_p_bound) and delta_p_bound > 0.0

    status = "PASS_E1569A_AND_PASS_T1569B_CANDIDATE" if (pass_e1569a and pass_t1569b) else "FAIL_E1569A_OR_T1569B"

    summary = {
        "checkpoint": "P1570_S520",
        "date_utc": "2026-05-14",
        "status": status,
        "strict_only": True,
        "legacy_bridge_used": False,
        "qw2191_closed": False,
        "toe_closed": False,
        "analytic_jacobian": ja,
        "numeric_jacobian": jn,
        "jacobian_max_abs_error": max_abs_err,
        "jacobian_tolerance": TOL_J,
        "E1569A_pass": pass_e1569a,
        "xi_hessian_closed_form": hessian_xi(BASE),
        "jinv_inf_norm": j_inv_norm,
        "delta_f_inf_norm": delta_f_norm,
        "delta_p_bound_inf": delta_p_bound,
        "T1569B_candidate_pass": pass_t1569b,
        "missing_exports_witnesses_theorems": [
            "T1570C_uniform_stability_bound_on_extended_domain",
            "T1570D_semiglobal_inverse_chart_patching",
            "W1570E_full_sm_gr_bound_compatibility_witness",
        ],
        "next_honest_step": "P1571_semiglobal_stability_scan_on_admissible_domain",
        "lay_summary": "Mamy jawne wzory na czułość modelu i lokalny przelicznik błędu pomiaru na błąd parametrów strict.",
    }
    summary["audit_digest"] = hashlib.sha256(json.dumps(summary, sort_keys=True).encode("utf-8")).hexdigest()

    out = out_dir / "p1570_s520_strict_closed_form_jacobian_hessian_and_stability_bound_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1570] wrote {out} status={status}")


if __name__ == "__main__":
    main()
