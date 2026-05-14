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

    good, buffer, bad = [], [], []
    norms = []
    for p in itertools.product(OMEGA_GRID, PHI_GRID, BETA_GRID, ETA_GRID):
        jinv = inv4(jacobian4(p))
        ninf = inf_norm(jinv)
        norms.append(ninf)
        entry = {"omega": p[0], "phi": p[1], "beta": p[2], "eta": p[3], "jinv_inf": ninf}
        if ninf <= GOOD_CUT:
            good.append(entry)
        elif ninf <= BUFFER_CUT:
            buffer.append(entry)
        else:
            bad.append(entry)

    b_good = max(x["jinv_inf"] for x in good) if good else float("inf")
    b_buffer = max(x["jinv_inf"] for x in buffer) if buffer else float("inf")
    b_bad = max(x["jinv_inf"] for x in bad) if bad else float("inf")

    cond1 = len(good) > 0 and len(buffer) > 0
    cond2 = len(bad) > 0
    cond3 = math.isfinite(b_good) and math.isfinite(b_buffer) and math.isfinite(b_bad) and (b_good < b_buffer < b_bad)

    status = "PASS_T1571A_CANDIDATE" if (cond1 and cond2 and cond3) else "FAIL_T1571A_CANDIDATE"

    summary = {
        "checkpoint": "P1572_S522",
        "date_utc": "2026-05-14",
        "status": status,
        "strict_only": True,
        "legacy_bridge_used": False,
        "qw2191_closed": False,
        "toe_closed": False,
        "grid_points": len(norms),
        "cuts": {"good": GOOD_CUT, "buffer": BUFFER_CUT},
        "chart_sizes": {"good": len(good), "buffer": len(buffer), "bad": len(bad)},
        "regional_bounds": {"B_good": b_good, "B_buffer": b_buffer, "B_bad": b_bad},
        "conditions": {
            "good_and_buffer_nonempty": cond1,
            "bad_chart_separated": cond2,
            "monotone_bounds": cond3,
        },
        "worst_bad_point": max(bad, key=lambda x: x["jinv_inf"]) if bad else None,
        "missing_exports_witnesses_theorems": [
            "T1572B_chart_transition_continuity_theorem",
            "W1572C_sm_gr_bundle_patching_compatibility_witness",
            "T1572D_global_inverse_error_theorem_for_eom_measurement_class",
        ],
        "next_honest_step": "P1573_formalize_chart_transition_maps_and_continuity",
        "lay_summary": "Domena została uczciwie podzielona na strefy stabilności; każda strefa ma własny, jawny limit błędu rekonstrukcji.",
    }
    summary["audit_digest"] = hashlib.sha256(json.dumps(summary, sort_keys=True).encode("utf-8")).hexdigest()

    out = out_dir / "p1572_s522_strict_semiglobal_chart_patching_and_regional_bound_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1572] wrote {out} status={status}")


if __name__ == "__main__":
    main()
