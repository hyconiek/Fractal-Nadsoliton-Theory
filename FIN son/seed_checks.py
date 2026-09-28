"""seed_checks.py -- reproduction lock + the probes that drove the architecture decisions.
Usage:  python3 seed_checks.py [--with-n6]      (N<=5 by default; N=6 costs ~2-3 min per model)
Writes seed_results.json next to this file.  Every 'repo_*' value is quoted from the repository's own reports
and is used ONLY as a reproduction target (never fitted)."""
import json
import math
import sys
import time

import numpy as np

from fin_core import (A7_default, G_FROZEN, A_from_lams, V_d4_saddle, V_localized, g_eq,
                      projector_countermodel, slow_modes, static_S)

R = {}
A7, L = A7_default()
Ap, lbar = projector_countermodel()
G = G_FROZEN

# 1. ---- reproduction lock (repo targets quoted from FIN_FULL_HANDOFF_327_CURRENT_20260928 / 300_326 reports)
V_loc = V_localized(G, A7)[0]
V_d4 = V_d4_saddle(G, A7)[0]
R["repro"] = {
    "L_k(k=1..6)": [float(x) for x in L[1:7]],
    "A7_first_row": [float(x) for x in A7[0]],
    "V_loc(G)": V_loc, "repo_V_loc": -0.80531946214230935605557878998,
    "V_d4(G)": V_d4, "repo_V_d4": -0.14310032501484691985100395454,
    "B4": V_d4 - V_loc, "repo_B4": 0.6622191371274624362045748354471982765,
    "g6=12/L6": 12 / L[6], "repo_g6": 5.1234275514,
    "spinodal_ladder_12/L_k(k=3..6)": [12 / float(L[k]) for k in (3, 4, 5, 6)],
    "g_eq": g_eq(A7, lo=2.0, hi=8.0, n=40), "G/g_eq": None, "G/g6": G / (12 / L[6]),
}
R["repro"]["G/g_eq"] = G / R["repro"]["g_eq"]
print("repro:", json.dumps({k: v for k, v in R["repro"].items() if not isinstance(v, list)}, indent=1), flush=True)

# 2. ---- rho_N (k=4 slow mode) vs repo 304 values
repo_rho = {3: 0.131439786195646, 4: 0.071692277481914, 5: 0.0401867559160965, 6: 0.0226112751871144}
Ns = [3, 4, 5] + ([6] if "--with-n6" in sys.argv else [])
R["rho_repro"] = {}
for N in Ns:
    rho, Rk, asym, fast = slow_modes(N, G, A7, "heatbath")
    R["rho_repro"][N] = {"rho": rho, "repo_rho": repo_rho[N], "rel_err": abs(rho - repo_rho[N]) / repo_rho[N], "DB_asym": asym}
    print("rho", N, R["rho_repro"][N], flush=True)

# 3. ---- cutoff table: coexistence gain for a projector with unit eigenvalue on harmonic set S
def Aproj(modes):
    return A_from_lams({k: 1.0 for k in modes})

R["cutoff_gEq_unit_eigenvalue"] = {}
for name, S in [("Potts k=1..6", [1, 2, 3, 4, 5, 6]), ("k>=2", [2, 3, 4, 5, 6]), ("k>=3 (FIN cutoff)", [3, 4, 5, 6]),
                ("k>=4", [4, 5, 6]), ("k>=5", [5, 6])]:
    R["cutoff_gEq_unit_eigenvalue"][name] = g_eq(Aproj(S), lo=3.0, hi=20.0, n=60)
    print("cutoff", name, R["cutoff_gEq_unit_eigenvalue"][name], flush=True)
R["cutoff_gEq_unit_eigenvalue"]["textbook Potts 2(q-1)ln(q-1)/(q-2)"] = 2 * 11 * math.log(11) / 10

# 4. ---- statics: projector vs A7
g_p = G * g_eq(Ap, lo=2.0, hi=8.0, n=40) / R["repro"]["g_eq"]
R["projector"] = {"lbar": lbar, "g_eq_proj": g_p / G * R["repro"]["g_eq"], "g_eq*lbar": g_eq(Ap, lo=2.0, hi=8.0, n=40) * lbar,
                  "g_eq_FIN*lbar": R["repro"]["g_eq"] * lbar,
                  "B4_proj_at_same_G/g_eq": V_d4_saddle(g_p, Ap)[0] - V_localized(g_p, Ap)[0]}
print("projector:", R["projector"], flush=True)

# 5. ---- R_k under three kinetics, FIN vs projector (same G/g_eq)
R["Rk"] = {}
for N in [n for n in Ns if n in (4, 5, 6)]:
    for kind in ("heatbath", "metropolis", "barker"):
        t = time.time()
        rf, Rf, asym, fast = slow_modes(N, G, A7, kind)
        rp, Rp, _, _ = slow_modes(N, g_p, Ap, kind)
        R["Rk"][f"N{N}|{kind}"] = {"rho_FIN": rf, "R_FIN": Rf, "R_proj": Rp, "DB_asym": asym, "fast_band_start": fast,
                                   "rel_diff_pct_FIN_vs_proj": {k: 100 * (Rf[k] - Rp[k]) / Rp[k] for k in Rf if k != 4}}
        print(f"Rk N={N} {kind}: R_FIN=" + " ".join(f"k{k}:{Rf[k]:.4f}" for k in sorted(Rf)) +
              "  FINvsProj%=" + " ".join(f"k{k}:{R['Rk'][f'N{N}|{kind}']['rel_diff_pct_FIN_vs_proj'][k]:+.1f}" for k in sorted(Rf) if k != 4) +
              f" [{time.time() - t:.0f}s]", flush=True)

# 6. ---- exact finite-N static Fourier variances (NOT Gaussian) FIN vs projector
R["static_S"] = {}
for N in (6, 8):
    for g in (2.0, 3.0):
        R["static_S"][f"N{N}|g{g}"] = {"FIN": static_S(N, g, A7), "proj_same_g": static_S(N, g, Ap),
                                       "gaussian_pred_FIN": {k: (1 / (1 - g * float(L[k]) / 12) if k >= 3 else 1.0) for k in range(1, 7)}}
        print("static", N, g, {k: round(v, 3) for k, v in R["static_S"][f"N{N}|g{g}"]["FIN"].items()}, flush=True)

json.dump(R, open("seed_results.json", "w"), indent=1, default=float)
print("written seed_results.json")
