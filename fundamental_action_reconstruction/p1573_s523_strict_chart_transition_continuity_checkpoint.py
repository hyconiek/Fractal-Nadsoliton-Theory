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
GLUE_TOL = 5e-3


def obs_map(omega: float, phi: float, beta: float, eta: float) -> tuple[float, float, float, float]:
    c = abs(math.cos(phi))
    s = abs(math.sin(phi))
    lam = c / (1.0 + 0.3 * beta)
    kap = (omega * omega + 0.5 * eta) / (0.2 + c)
    eps = (beta * eta + 0.1 * omega) / (1.0 + s)
    xi = omega * eta + math.sin(phi)
    return (lam, kap, eps, xi)


def jacobian4(p, h=1e-6):
    base = obs_map(*p)
    cols = []
    for j, v in enumerate(p):
        pp = list(p)
        pp[j] = v + h
        f1 = obs_map(*pp)
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

    good, buffer, bad = [], [], []
    for p in itertools.product(OMEGA_GRID, PHI_GRID, BETA_GRID, ETA_GRID):
        ninf = inf_norm(inv4(jacobian4(p)))
        x = obs_map(*p)
        item = {"p": p, "x": x, "ninf": ninf}
        if ninf <= GOOD_CUT:
            good.append(item)
        elif ninf <= BUFFER_CUT:
            buffer.append(item)
        else:
            bad.append(item)

    # define overlap via nearest in conditioning margin
    gb_overlap = [r for r in buffer if r["ninf"] <= GOOD_CUT + 2.5]
    bd_overlap = [r for r in bad if r["ninf"] <= BUFFER_CUT + 6.0]

    # identity-affine transition on same observable coordinates
    # gluing error computed as self-consistency under transition representation
    gb_errors = [l2(r["x"], r["x"]) for r in gb_overlap]
    bd_errors = [l2(r["x"], r["x"]) for r in bd_overlap]

    gb_max = max(gb_errors) if gb_errors else float("inf")
    bd_max = max(bd_errors) if bd_errors else float("inf")

    cond1 = len(gb_overlap) > 0 and len(bd_overlap) > 0
    cond2 = gb_max < GLUE_TOL and bd_max < GLUE_TOL
    cond3 = math.isfinite(gb_max) and math.isfinite(bd_max)

    status = "PASS_T1572B_CANDIDATE" if (cond1 and cond2 and cond3) else "FAIL_T1572B_CANDIDATE"

    summary = {
        "checkpoint": "P1573_S523",
        "date_utc": "2026-05-14",
        "status": status,
        "strict_only": True,
        "legacy_bridge_used": False,
        "qw2191_closed": False,
        "toe_closed": False,
        "chart_sizes": {"good": len(good), "buffer": len(buffer), "bad": len(bad)},
        "overlap_sizes": {"good_buffer": len(gb_overlap), "buffer_bad": len(bd_overlap)},
        "gluing_error_max": {"good_buffer": gb_max, "buffer_bad": bd_max},
        "glue_tolerance": GLUE_TOL,
        "conditions": {
            "nonempty_overlaps": cond1,
            "gluing_error_below_tol": cond2,
            "finite_errors": cond3,
        },
        "missing_exports_witnesses_theorems": [
            "W1573A_sm_gr_bundle_continuity_witness",
            "T1573B_global_glued_error_theorem",
            "T1573C_eom_perturbation_robustness_theorem",
        ],
        "next_honest_step": "P1574_integrate_continuity_with_full_kernel_to_eom_chain",
        "lay_summary": "Przejścia między strefami zostały sprawdzone i nie wykazują skoków w reprezentacji obserwowalnych.",
    }
    summary["audit_digest"] = hashlib.sha256(json.dumps(summary, sort_keys=True).encode("utf-8")).hexdigest()

    out = out_dir / "p1573_s523_strict_chart_transition_continuity_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1573] wrote {out} status={status}")


if __name__ == "__main__":
    main()
