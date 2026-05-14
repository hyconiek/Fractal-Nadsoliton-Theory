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
B_GLOBAL_THRESHOLD = 40.0


def coeffs(o, ph, b, e):
    lam = abs(math.cos(ph)) / (1.0 + 0.3 * b)
    kap = (o * o + 0.5 * e) / (0.2 + abs(math.cos(ph)))
    eps = (b * e + 0.1 * o) / (1.0 + abs(math.sin(ph)))
    return (lam, kap, eps)


def chain_gain(c):
    lam, kap, eps = c
    return math.sqrt((2 * lam - eps) ** 2 + (2 * kap + 0.5 * eps) ** 2 + (eps - lam) ** 2)


def jac4(p, h=1e-6):
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


def main() -> None:
    root = Path(__file__).resolve().parent
    out_dir = root / "generated"
    out_dir.mkdir(parents=True, exist_ok=True)

    overlap = []
    for p in itertools.product(OMEGA_GRID, PHI_GRID, BETA_GRID, ETA_GRID):
        ninf = inf_norm(inv4(jac4(p)))
        if GOOD_CUT < ninf <= BUFFER_CUT + 6.0:
            c = coeffs(*p)
            gain = chain_gain(c)
            b_local = ninf * gain
            overlap.append({"p": p, "jinv_inf": ninf, "gain": gain, "B_local": b_local})

    overlap.sort(key=lambda r: r["p"])
    b_global = max(r["B_local"] for r in overlap) if overlap else float("inf")
    top = sorted(overlap, key=lambda r: r["B_local"], reverse=True)[:5]

    status = "PASS_T1574B_CANDIDATE" if math.isfinite(b_global) and b_global < B_GLOBAL_THRESHOLD else "FAIL_T1574B_CANDIDATE"

    summary = {
        "checkpoint": "P1575_S525",
        "date_utc": "2026-05-14",
        "status": status,
        "strict_only": True,
        "legacy_bridge_used": False,
        "qw2191_closed": False,
        "toe_closed": False,
        "overlap_points": len(overlap),
        "B_global": b_global,
        "B_global_threshold": B_GLOBAL_THRESHOLD,
        "top_dominant_points": [
            {"omega": r["p"][0], "phi": r["p"][1], "beta": r["p"][2], "eta": r["p"][3], "B_local": r["B_local"]}
            for r in top
        ],
        "missing_exports_witnesses_theorems": [
            "T1575A_eom_perturbation_robustness_theorem",
            "W1575B_toe_target_alignment_witness",
            "T1575C_final_global_gluing_theorem",
        ],
        "next_honest_step": "P1576_stress_test_eom_perturbations_and_formalize_T1575A",
        "lay_summary": "Wyznaczono globalny licznik wzmacniania błędu i punkty, które najbardziej go napędzają.",
    }
    summary["audit_digest"] = hashlib.sha256(json.dumps(summary, sort_keys=True).encode("utf-8")).hexdigest()

    out = out_dir / "p1575_s525_strict_global_error_propagation_bound_summary.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1575] wrote {out} status={status}")


if __name__ == "__main__":
    main()
