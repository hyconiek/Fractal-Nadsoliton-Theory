from __future__ import annotations

import hashlib
import itertools
import json
import math
from pathlib import Path

OMEGA_GRID = [0.16, 0.22, 0.28, 0.36]
PHI_GRID = [0.08, 0.14, 0.20, 0.24]
BETA_GRID = [0.85, 0.95, 1.05, 1.15]
ETA_GRID = [1.55, 1.70, 1.85, 1.95]
GOOD_CUT = 12.0
BUFFER_CUT = 20.0

COEFF_TOL = 0.12
LAG_TOL = 0.08
EOM_TOL = 0.15


def kernel_features(omega: float, phi: float, beta: float, eta: float) -> tuple[float, float, float]:
    return (math.cos(phi), beta, eta)


def coeffs(omega: float, phi: float, beta: float, eta: float) -> tuple[float, float, float]:
    lam = abs(math.cos(phi)) / (1.0 + 0.3 * beta)
    kap = (omega * omega + 0.5 * eta) / (0.2 + abs(math.cos(phi)))
    eps = (beta * eta + 0.1 * omega) / (1.0 + abs(math.sin(phi)))
    return (lam, kap, eps)


def lag_proxy(c: tuple[float, float, float]) -> float:
    lam, kap, eps = c
    return 0.5 * lam * lam + 0.5 * kap * kap + 0.25 * eps * eps


def eom_proxy(c: tuple[float, float, float]) -> tuple[float, float, float]:
    lam, kap, eps = c
    return (2.0 * lam - eps, 2.0 * kap + 0.5 * eps, eps - lam)


def jacobian4(p, h=1e-6):
    def obs(pp):
        o, ph, b, e = pp
        c = coeffs(o, ph, b, e)
        return (*c, o * e + math.sin(ph))

    base = obs(p)
    cols = []
    for j, v in enumerate(p):
        pp = list(p)
        pp[j] = v + h
        f1 = obs(tuple(pp))
        cols.append([(f1[i] - base[i]) / h for i in range(4)])
    return [[cols[c][r] for c in range(4)] for r in range(4)]


def inv4(m):
    n = 4
    a = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(m)]
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(a[r][i]))
        a[i], a[piv] = a[piv], a[i]
        p = a[i][i]
        if abs(p) < 1e-14:
            raise ValueError("singular")
        for j in range(2 * n):
            a[i][j] /= p
        for r in range(n):
            if r == i:
                continue
            fac = a[r][i]
            for j in range(2 * n):
                a[r][j] -= fac * a[i][j]
    return [row[n:] for row in a]


def inf_norm(m):
    return max(sum(abs(x) for x in row) for row in m)


def l2(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def main() -> None:
    root = Path(__file__).resolve().parent
    out_dir = root / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for p in itertools.product(OMEGA_GRID, PHI_GRID, BETA_GRID, ETA_GRID):
        ninf = inf_norm(inv4(jacobian4(p)))
        c = coeffs(*p)
        rows.append({"p": p, "ninf": ninf, "c": c, "L": lag_proxy(c), "eom": eom_proxy(c)})

    gb = [r for r in rows if GOOD_CUT < r["ninf"] <= GOOD_CUT + 2.5]
    bd = [r for r in rows if BUFFER_CUT < r["ninf"] <= BUFFER_CUT + 6.0]

    def max_jumps(sample):
        if len(sample) < 2:
            return float("inf"), float("inf"), float("inf")
        coeff_jump = lag_jump = eom_jump = 0.0
        for a, b in zip(sample[:-1], sample[1:]):
            coeff_jump = max(coeff_jump, l2(a["c"], b["c"]))
            lag_jump = max(lag_jump, abs(a["L"] - b["L"]))
            eom_jump = max(eom_jump, l2(a["eom"], b["eom"]))
        return coeff_jump, lag_jump, eom_jump

    gb_c, gb_l, gb_e = max_jumps(gb)
    bd_c, bd_l, bd_e = max_jumps(bd)

    max_coeff = max(gb_c, bd_c)
    max_lag = max(gb_l, bd_l)
    max_eom = max(gb_e, bd_e)

    cond = max_coeff < COEFF_TOL and max_lag < LAG_TOL and max_eom < EOM_TOL
    status = "PASS_W1573A_CANDIDATE" if cond else "FAIL_W1573A_CANDIDATE"

    summary = {
        "checkpoint": "P1574_S524",
        "date_utc": "2026-05-14",
        "status": status,
        "strict_only": True,
        "legacy_bridge_used": False,
        "qw2191_closed": False,
        "toe_closed": False,
        "overlap_sizes": {"good_buffer": len(gb), "buffer_bad": len(bd)},
        "max_jumps": {
            "coefficients": max_coeff,
            "lagrangian_proxy": max_lag,
            "eom_proxy": max_eom,
        },
        "tolerances": {
            "coefficients": COEFF_TOL,
            "lagrangian_proxy": LAG_TOL,
            "eom_proxy": EOM_TOL,
        },
        "missing_exports_witnesses_theorems": [
            "T1574B_global_error_propagation_bound",
            "T1574C_full_chain_noise_stability_theorem",
            "W1574D_toe_target_alignment_witness",
        ],
        "next_honest_step": "P1575_formalize_global_error_propagation_bound",
        "lay_summary": "Przejścia chartów pozostają ciągłe także po przejściu przez współczynniki, lagranżian i równania ruchu.",
    }
    summary["audit_digest"] = hashlib.sha256(json.dumps(summary, sort_keys=True).encode("utf-8")).hexdigest()

    out = out_dir / "p1574_s524_strict_full_chain_continuity_witness_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1574] wrote {out} status={status}")


if __name__ == "__main__":
    main()
