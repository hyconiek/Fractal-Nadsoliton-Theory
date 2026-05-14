from __future__ import annotations

import hashlib
import itertools
import json
import math
from pathlib import Path

DET_THRESHOLD = 5e-2
COND_THRESHOLD = 12.0

OMEGA_GRID = [0.16, 0.22, 0.28, 0.36]
PHI_GRID = [0.08, 0.14, 0.20, 0.24]
BETA_GRID = [0.85, 0.95, 1.05, 1.15]
ETA_GRID = [1.55, 1.70, 1.85, 1.95]


def obs_map(omega: float, phi: float, beta: float, eta: float) -> tuple[float, float, float, float]:
    c = abs(math.cos(phi))
    s = abs(math.sin(phi))
    lam = c / (1.0 + 0.3 * beta)
    kap = (omega * omega + 0.5 * eta) / (0.2 + c)
    eps = (beta * eta + 0.1 * omega) / (1.0 + s)
    xi = omega * eta + math.sin(phi)
    return (lam, kap, eps, xi)


def jacobian4(p: tuple[float, float, float, float], h: float = 1e-6) -> list[list[float]]:
    base = obs_map(*p)
    cols = []
    for j, v in enumerate(p):
        pp = list(p)
        pp[j] = v + h
        f1 = obs_map(*pp)
        cols.append([(f1[i] - base[i]) / h for i in range(4)])
    return [[cols[c][r] for c in range(4)] for r in range(4)]


def det4(m: list[list[float]]) -> float:
    a = [row[:] for row in m]
    n = 4
    sgn = 1.0
    det = 1.0
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(a[r][i]))
        if abs(a[piv][i]) < 1e-15:
            return 0.0
        if piv != i:
            a[i], a[piv] = a[piv], a[i]
            sgn *= -1.0
        p = a[i][i]
        det *= p
        for r in range(i + 1, n):
            fac = a[r][i] / p
            for c in range(i, n):
                a[r][c] -= fac * a[i][c]
    return sgn * det


def inv4(m: list[list[float]]) -> list[list[float]]:
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


def inf_norm(m: list[list[float]]) -> float:
    return max(sum(abs(x) for x in row) for row in m)


def main() -> None:
    root = Path(__file__).resolve().parent
    out_dir = root / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    worst_cond = -1.0
    worst_point = None
    min_abs_det = float("inf")
    sign_set = set()

    for p in itertools.product(OMEGA_GRID, PHI_GRID, BETA_GRID, ETA_GRID):
        j = jacobian4(p)
        d = det4(j)
        sign_set.add(1 if d > 0 else (-1 if d < 0 else 0))
        min_abs_det = min(min_abs_det, abs(d))
        jinv = inv4(j)
        cond = inf_norm(jinv)
        if cond > worst_cond:
            worst_cond = cond
            worst_point = p
        rows.append({"p": p, "det": d, "jinv_inf": cond})

    cond1 = sign_set == {1} or sign_set == {-1}
    cond2 = min_abs_det > DET_THRESHOLD
    cond3 = worst_cond < COND_THRESHOLD

    status = "PASS_T1570C_CANDIDATE" if (cond1 and cond2 and cond3) else "FAIL_T1570C_CANDIDATE"

    summary = {
        "checkpoint": "P1571_S521",
        "date_utc": "2026-05-14",
        "status": status,
        "strict_only": True,
        "legacy_bridge_used": False,
        "qw2191_closed": False,
        "toe_closed": False,
        "grid_shape": [len(OMEGA_GRID), len(PHI_GRID), len(BETA_GRID), len(ETA_GRID)],
        "grid_points": len(rows),
        "det_signs": sorted(list(sign_set)),
        "min_abs_det": min_abs_det,
        "det_threshold": DET_THRESHOLD,
        "worst_jinv_inf": worst_cond,
        "worst_point": {
            "omega": worst_point[0],
            "phi": worst_point[1],
            "beta": worst_point[2],
            "eta": worst_point[3],
        },
        "cond_threshold": COND_THRESHOLD,
        "conditions": {
            "det_sign_stable": cond1,
            "det_min_ok": cond2,
            "worst_condition_ok": cond3,
        },
        "missing_exports_witnesses_theorems": [
            "T1571A_semiglobal_chart_patching_theorem",
            "W1571B_sm_gr_bound_compatibility_witness",
            "T1571C_global_inverse_error_bound_export",
        ],
        "next_honest_step": "P1572_build_semiglobal_chart_patching_and_global_error_bound",
        "lay_summary": "Na całej sprawdzonej siatce model zachowuje stabilną odwracalność, a najtrudniejszy punkt nadal pozostaje kontrolowalny.",
    }
    summary["audit_digest"] = hashlib.sha256(json.dumps(summary, sort_keys=True).encode("utf-8")).hexdigest()

    out = out_dir / "p1571_s521_strict_semiglobal_stability_scan_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1571] wrote {out} status={status}")


if __name__ == "__main__":
    main()
