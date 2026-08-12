#!/usr/bin/env python3
"""Live and serialized checks for FIN ST263--ST277."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from fin_st01_st15_research import N, strict_operator
from fin_st263_st277_research import (
    PACKETS, build_local_refinement, d12_permutations, real_fourier_modes,
    selector_float, selector_geometry,
)

ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / "FIN_ST263_ST277_Results.json").read_text())
checks = []


def check(name, condition):
    checks.append((name, bool(condition)))


for k in range(263, 278):
    p = PACKETS[k]
    check(f"ST{k} packet exists", p.exists())
    check(f"ST{k} packet hash", hashlib.sha256(p.read_bytes()).hexdigest() == DATA[f"ST{k}"]["packet_sha256"])

w, a, _ = strict_operator()
ev = np.linalg.eigvalsh(a)
check("ST263 lowest multiplicity", np.sum(np.isclose(ev, ev[1], atol=1e-11)) == 2)
check("ST263 degree-12 invariant", DATA["ST263"]["first_anisotropic_invariant"] == "Re(z^12)")
check("ST263 twelve branches", DATA["ST263"]["number_of_phase_branches"] == 12)
check("ST263 ensemble uniform", DATA["ST263"]["ensemble_distance_from_uniform"] < 1e-14)
check("ST263 states positive", DATA["ST263"]["branch_minimum_probability"] > 0)

mons = {(x["p"], x["q"], x["degree"]) for x in DATA["ST264"]["equivariant_monomials"]}
expected = {(q + 1, q, 2*q + 1) for q in range(6)} | {(0, 11, 11)}
check("ST264 module basis", mons == expected)
check("ST264 dimension", DATA["ST264"]["real_module_dimension"] == 7)
check("ST264 covariance", DATA["ST264"]["maximum_numeric_covariance_error"] < 1e-15)

check("ST265 selected pair", DATA["ST265"]["selected_seed_pair"] == [20, 23])
check("ST265 44 boxes", sum(DATA["ST265"]["new_steps_per_seed"].values()) == 44)
check("ST265 separation", DATA["ST265"]["minimum_new_pair_clearance"] > 0)
check("ST265 interval rank", DATA["ST265"]["minimum_interval_full_row_rank_lower_bound"] > 0)
for rows in DATA["ST265"]["paths"].values():
    check("ST265 all Krawczyk inclusions", all(r["certificate"]["included"] for r in rows))

st266 = DATA["ST266"]
check("ST266 optimum", st266["entanglement_fidelity_optimum"] == "3/8")
check("ST266 zero gap", abs(st266["primal_discrimination_sum"] - st266["dual_trace"]) < 1e-14)
check("ST266 PSD slack0", min(st266["dual_slack_0_eigenvalues"]) > -1e-14)
check("ST266 PSD slack1", min(st266["dual_slack_1_eigenvalues"]) > -1e-14)
check("ST266 nonunital", st266["unitality_defect_Frobenius"] > .1)
check("ST266 noncommuting outputs", st266["prepared_state_commutator_norm"] > .1)

st267 = DATA["ST267"]
z = np.array(st267["KKT_center"])
_, g, h, afinal, cfinal = selector_float(z[:10])
_, _, _, b, bc, _, dtk = selector_geometry()
F = np.r_[g + z[10]*afinal - z[11]*b, cfinal + afinal@z[:10] - .6, b@z[:10] + bc]
check("ST267 live KKT residual", np.linalg.norm(F, np.inf) < 1e-12)
check("ST267 Krawczyk", st267["minimum_Krawczyk_margin"] > 0)
check("ST267 positive multiplier", st267["active_multiplier"] > 0)
check("ST267 Hessian positive", np.linalg.eigvalsh(h)[0] > .4)
check("ST267 inactive slacks", st267["minimum_inactive_slew_margin"] > .03)
check("ST267 endpoint", abs(cfinal + afinal@z[:10] - .6) < 1e-12)
check("ST267 active final slew", abs((z[9]-z[8])/dtk + .3) < 1e-12)

st268 = DATA["ST268"]
check("ST268 intertwining", st268["intertwining_error"] < 1e-12)
check("ST268 metric naturality", st268["metric_naturality_error"] < 1e-12)
check("ST268 mean-zero", st268["mean_zero_transport_error"] < 1e-12)

st269 = DATA["ST269"]
check("ST269 full count", st269["last_complete_boxes"] == 64000)
check("ST269 no failure at .00070", st269["last_complete_failures"] == 0)
check("ST269 positive full margin", st269["last_complete_minimum_margin"] > 0)
check("ST269 failures at .00075", st269["next_attempt_failed_corner_cells"] == 5)

modes = real_fourier_modes()
check("ST270 modes orthonormal", np.linalg.norm(modes.T@modes-np.eye(6)) < 1e-12)
check("ST270 coverage", all(abs(r["simultaneous_confidence"]-.95) < 1e-15 for r in DATA["ST270"]["rows"]))
check("ST270 design factor", all(abs(r["exact_covariance_design_factor"]-1.95) < 1e-14 for r in DATA["ST270"]["rows"]))
check("ST270 finite widths", all(np.isfinite(r["maximum_halfwidth"]) and r["maximum_halfwidth"] > 0 for r in DATA["ST270"]["rows"]))

st271 = DATA["ST271"]
check("ST271 exact family", max(r["involution_error"] for r in st271["sampled_family"]) < 1e-14)
check("ST271 anticommutes", max(r["anticommutator_error"] for r in st271["sampled_family"]) < 1e-14)
check("ST271 constant objective", np.ptp([r["objective"] for r in st271["sampled_family"]]) < 1e-14)

st272 = DATA["ST272"]
widths = [r["spectral_radius_interval"][1]-r["spectral_radius_interval"][0] for r in st272["adaptive_history"]]
check("ST272 adaptive widths decrease", widths[2] < widths[1] < widths[0])
check("ST272 fixed s rigorous", st272["certified_fixed_s"] == "1/2")
check("ST272 Chernoff not promoted", st272["optimization_status"] == "strong numerical evidence only")
check("ST272 combined improves", st272["combined_width"] < .10885)

check("ST273 controlled unitary", DATA["ST273"]["controlled_SWAP_unitarity_error"] < 1e-14)
check("ST273 involution", DATA["ST273"]["controlled_SWAP_involution_error"] < 1e-14)
check("ST273 entropy nonnegative", all(r["entropy_nats"] >= 0 for r in DATA["ST273"]["degenerate_controller_rows"]))

for split in (.5,):
    at = build_local_refinement(w, split, .7)
    plus = np.array([1,1])/math.sqrt(2); R=np.kron(np.eye(N),plus[:,None])
    check("ST274 live intertwining", np.linalg.norm(at@R-R@a) < 1e-12)
check("ST274 dimension reduction", DATA["ST274"]["dimension_before"] == 7 and DATA["ST274"]["dimension_after"] == 1)
check("ST274 relabel natural", max(r["single_fiber_relabeling_error"] for r in DATA["ST274"]["rows"]) < 1e-12)

st275 = DATA["ST275"]
check("ST275 bound", st275["maximum_bound_violation"] <= 0)
check("ST275 trace obstruction", "diverges" in st275["full_log_Haar_heat_trace_status"])
check("ST275 positive cutoff curve", min(st275["finite_cutoff_curve"]) >= 0)
check("ST275 below plateau", max(st275["finite_cutoff_curve"]) <= st275["plateau"] + 1e-14)

check("ST276 heat contrast", DATA["ST276"]["recovered_heat_contrast_error"] < 1e-12)
check("ST276 principal log", DATA["ST276"]["principal_log_B_reconstruction_error"] < 1e-12)
check("ST276 coarse heat", DATA["ST276"]["coarse_heat_error"] < 1e-12)
check("ST276 record fields", len(DATA["ST276"]["record_schema"]) == 14)
check("ST277 blocked", DATA["ST277"]["status"] == "blocked_no_external_record")
check("15 recommendations", len(DATA["recommended_next_programs"]) == 15)
check("epistemic no ToE", "ToE closure" in DATA["epistemic_boundary"])

failed = [name for name, ok in checks if not ok]
if failed:
    raise SystemExit("FAILED: " + "; ".join(failed))
print(f"{len(checks)}/{len(checks)} tests passed")
