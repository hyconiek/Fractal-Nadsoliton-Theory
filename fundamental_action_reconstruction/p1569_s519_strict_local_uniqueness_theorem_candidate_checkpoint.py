from __future__ import annotations

import hashlib
import itertools
import json
import math
from pathlib import Path

BASE = (0.18575, 0.16250, 1.0, 1.8)
H = 1e-6
RADIUS = 0.01
GRID = (-1.0, 0.0, 1.0)
DET_THRESHOLD = 1e-3


def obs_map(omega: float, phi: float, beta: float, eta: float) -> tuple[float, float, float, float]:
    lam = abs(math.cos(phi)) / (1.0 + 0.3 * beta)
    kap = (omega * omega + 0.5 * eta) / (0.2 + abs(math.cos(phi)))
    eps = (beta * eta + 0.1 * omega) / (1.0 + abs(math.sin(phi)))
    xi = omega * eta + math.sin(phi)
    return (lam, kap, eps, xi)


def det4(m: list[list[float]]) -> float:
    # Gaussian elimination determinant
    a = [row[:] for row in m]
    n = 4
    sign = 1.0
    det = 1.0
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(a[r][i]))
        if abs(a[piv][i]) < 1e-15:
            return 0.0
        if piv != i:
            a[i], a[piv] = a[piv], a[i]
            sign *= -1.0
        p = a[i][i]
        det *= p
        for r in range(i + 1, n):
            fac = a[r][i] / p
            for c in range(i, n):
                a[r][c] -= fac * a[i][c]
    return sign * det


def jacobian4(params: tuple[float, float, float, float], h: float = H) -> list[list[float]]:
    base = obs_map(*params)
    cols = []
    for j, p in enumerate(params):
        pert = list(params)
        pert[j] = p + h
        f1 = obs_map(*pert)
        cols.append([(f1[i] - base[i]) / h for i in range(4)])
    return [[cols[c][r] for c in range(4)] for r in range(4)]


def main() -> None:
    root = Path(__file__).resolve().parent
    out_dir = root / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    j0 = jacobian4(BASE)
    det0 = det4(j0)

    local_dets: list[float] = []
    for shifts in itertools.product(GRID, repeat=4):
        p = tuple(BASE[i] + RADIUS * shifts[i] for i in range(4))
        d = det4(jacobian4(p))
        local_dets.append(d)

    signs = {1 if d > 0 else (-1 if d < 0 else 0) for d in local_dets}
    min_abs_det = min(abs(d) for d in local_dets)

    cond1 = abs(det0) > DET_THRESHOLD
    cond2 = signs == {1} or signs == {-1}
    cond3 = min_abs_det > DET_THRESHOLD

    status = "PASS_T1568D_CANDIDATE" if (cond1 and cond2 and cond3) else "FAIL_T1568D_CANDIDATE"

    summary = {
        "checkpoint": "P1569_S519",
        "date_utc": "2026-05-14",
        "status": status,
        "strict_only": True,
        "legacy_bridge_used": False,
        "qw2191_closed": False,
        "toe_closed": False,
        "base_params": {"omega": BASE[0], "phi": BASE[1], "beta": BASE[2], "eta": BASE[3]},
        "det_jacobian_base": det0,
        "det_threshold": DET_THRESHOLD,
        "local_radius": RADIUS,
        "local_grid_points": len(local_dets),
        "local_det_signs": sorted(list(signs)),
        "local_det_min_abs": min_abs_det,
        "conditions": {
            "det_base_nonzero": cond1,
            "sign_stability_local": cond2,
            "min_abs_det_local": cond3,
        },
        "missing_exports_witnesses_theorems": [
            "E1569A_closed_form_jacobian_hessian_export",
            "T1569B_inverse_error_stability_envelope",
            "T1569C_global_chart_patching_theorem",
        ],
        "next_honest_step": "P1570_export_closed_form_jacobian_hessian_and_stability_bound",
        "lay_summary": "W pobliżu punktu pracy 4 sygnały niosą wystarczającą informację, by odróżniać parametry strict jednoznacznie lokalnie.",
    }
    summary["audit_digest"] = hashlib.sha256(json.dumps(summary, sort_keys=True).encode("utf-8")).hexdigest()

    out = out_dir / "p1569_s519_strict_local_uniqueness_theorem_candidate_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1569] wrote {out} status={status}")


if __name__ == "__main__":
    main()
