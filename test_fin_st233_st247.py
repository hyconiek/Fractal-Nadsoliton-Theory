#!/usr/bin/env python3
"""Lightweight live/serialized validation for FIN programs ST233--ST247."""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from fin_st01_st15_research import N, strict_operator
from fin_st233_st247_research import optimal_equilibrium


ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / "FIN_ST233_ST247_Results.json").read_text(encoding="utf-8"))
COUNT = 0


def check(condition: bool, message: str) -> None:
    global COUNT
    if not condition:
        raise AssertionError(message)
    COUNT += 1


def main() -> None:
    check(DATA["metadata"]["programs"] == "ST233-ST247", "metadata range")
    check(all(f"ST{k}" in DATA for k in range(233, 248)), "complete program range")
    for k in range(233, 248):
        row = DATA[f"ST{k}"]
        path = ROOT / row["packet_file"]
        check(path.exists(), f"ST{k} packet exists")
        check(hashlib.sha256(path.read_bytes()).hexdigest() == row["packet_sha256"], f"ST{k} packet hash")

    _, a, _ = strict_operator()
    rho = expm(-0.7 * a)[:, 0]
    b = np.diag(rho)
    check(abs(rho.sum() - 1) < 2e-15 and rho.min() > 0, "ST233 live heat state")
    check(np.linalg.norm(a @ b - b @ a) > 0.25, "ST233 live noncommutation")
    check(DATA["ST233"]["uniform_state_commutator_norm"] == 0.0, "ST233 symmetric control")
    check(DATA["ST233"]["D12_covariance_error"] == 0.0, "ST233 covariance")

    st234 = DATA["ST234"]
    ratios = sorted(z["interval"] for z in st234["likelihood_ratio_intervals"])
    check(st234["all_heat_entries_strictly_positive"], "ST234 positivity")
    check(all(ratios[i][1] < ratios[i + 1][0] for i in range(11)), "ST234 disjoint intervals")
    check(st234["minimum_adjacent_interval_gap"] > 0.024, "ST234 separation margin")

    stab = {z["condition"]: z["stabilizer_order"] for z in DATA["ST235"]["rows"]}
    check(stab == {"uniform": 24, "localized delta_0": 2, "even strict heat state": 2, "chiral localized state": 1}, "ST235 stabilizers")
    check(DATA["ST235"]["rows"][-1]["vertex_orbits_under_stabilizer"] == 12, "ST235 chiral orbits")

    st236 = DATA["ST236"]
    check(st236["signed_seed_count"] == 60 and st236["certified_steps"] == 120, "ST236 bounded campaign")
    check(st236["finite_graph_component_count"] == 30 and len(st236["exact_uniform_fold_edges"]) == 30, "ST236 finite components")
    check(not st236["nonfold_box_intersections"] and st236["minimum_nonfold_box_clearance"] > 0.006, "ST236 separation")
    check(all(r["certificate"]["included"] for r in st236["rows"]), "ST236 inclusions")

    st237 = DATA["ST237"]
    check(sum(Fraction(x) for x in st237["exact_Pauli_probabilities"]) == 1, "ST237 probabilities")
    check(Fraction(st237["exact_optimum"]) == Fraction(189, 500), "ST237 exact optimum")
    check(st237["Bloch_determinant"] < 0 and abs(st237["primal_dual_gap"]) < 1e-15, "ST237 orientation and duality")
    check(min(st237["dual_slack_eigenvalues"]) > -1e-14, "ST237 dual PSD")

    st238 = DATA["ST238"]
    check(st238["time_interval_at_C_lower"][0] > 4 > st238["time_interval_at_C_upper"][1], "ST238 time bracket")
    check(st238["certified_C_bracket"][0] < st238["floating_C"] < st238["certified_C_bracket"][1], "ST238 C isolation")
    check(0.317 < st238["certified_entropy_production_bracket"][0] < st238["certified_entropy_production_bracket"][1] < 0.320, "ST238 cost")
    check(optimal_equilibrium(0.2, st238["floating_C"]) > 0.2, "ST238 live drift")

    st239 = DATA["ST239"]
    check(st239["nonlinear_chart_I12"] > st239["linear_chart_I12"], "ST239 chart obstruction")
    p = st239["parameters"]
    expected_shift = 6 * p["a"] * p["c"] ** 2 / p["a"] ** 6
    check(abs(st239["chart_induced_shift"] - expected_shift) < 1e-15, "ST239 live coefficient")

    st240 = DATA["ST240"]
    check(st240["coarse_cells_scanned"] == st240["failed_parents_subdivided"] == 64 and st240["coarse_passes"] == 0, "ST240 declared sample")
    check(st240["child_boxes"] == st240["certified_children"] == 512, "ST240 child resolution")
    check(st240["minimum_child_margin"] > 0.0017, "ST240 positive child margin")

    risks = DATA["ST241"]["rows"]
    check(abs(DATA["ST241"]["design_effect"] - 1.95) < 1e-14, "ST241 design effect")
    check(all(abs(r["cluster_trace_risk_constant"] / r["iid_trace_risk_constant"] - 1.95) < 1e-13 for r in risks), "ST241 cluster inflation")
    check(all(risks[i + 1]["iid_trace_risk_constant"] > risks[i]["iid_trace_risk_constant"] for i in range(len(risks) - 1)), "ST241 depth monotonicity")

    st242 = DATA["ST242"]
    check(all(r["minimum_singular_value_target_block"] > 0.4 for r in st242["rows"]), "ST242 invertible odd blocks")
    check(max(r["unitary_error"] for r in st242["rows"]) < 2e-15, "ST242 polar unitaries")
    check(max(r["anticommutator_error"] for r in st242["rows"]) == 0.0, "ST242 Clifford constraint")

    st243 = DATA["ST243"]
    check(st243["derivative_at_bracket"][0] < 0 < st243["derivative_at_bracket"][1], "ST243 numerical sign change")
    check(st243["optimized_block_coefficient_interval"][1] < st243["s_half_block_coefficient"][0], "ST243 improvement over s=1/2")
    check(st243["status"].startswith("strong_numerical"), "ST243 epistemic status")
    check(st243["optimized_asymptotic_Chernoff_rate_bracket"]["lower"] < st243["optimized_asymptotic_Chernoff_rate_bracket"]["upper"], "ST243 rate order")

    st244 = DATA["ST244"]
    check(all(r["swap_step_total_entropy_change"] == 0 for r in st244["rows"]), "ST244 SWAP entropy")
    check(abs(st244["rows"][1]["degenerate_H_blank_restoration_min_work_beta1"] - math.log(16)) < 1e-14, "ST244 Landauer recycling")

    st245 = DATA["ST245"]
    check([r["real_symmetric_moduli_dimension"] for r in st245["moduli_rows"]] == [78, 300, 666, 1176], "ST245 moduli dimensions")
    check(max(r["intertwining_error"] for r in st245["three_fiber_replays"]) < 4e-15, "ST245 intertwining replay")
    check(max(r["distance_replay_error_t1"] for r in st245["three_fiber_replays"]) < 1e-14, "ST245 metric replay")

    st246 = DATA["ST246"]
    check(len(st246["rows"]) == 7 and all(len(r["constant_curve"]) == 220 for r in st246["rows"]), "ST246 atlas dimensions")
    check(st246["rows"][-1]["constant_mu_peak"] > st246["rows"][0]["constant_mu_peak"], "ST246 refinement effect")
    check("t->t/c" in st246["scale_orbit"], "ST246 scale obstruction")

    check(DATA["ST247"]["status"] == "blocked_no_external_record" and not DATA["ST247"]["local_search_performed"], "ST247 evidence stop")
    check(len(DATA["recommended_next_programs"]) == 15 and DATA["recommended_next_programs"][0]["id"] == "ST248", "next recommendations")
    boundary = DATA["epistemic_boundary"]
    check("QW-2191" in boundary and "legacy-to-strict" in boundary and "ToE closure" in boundary, "epistemic guardrail")
    check("noncommuting second operator" in DATA["central_verdict"] and "nonuniqueness cone" in DATA["central_verdict"], "central verdict")

    print(f"{COUNT}/{COUNT} tests passed")


if __name__ == "__main__":
    main()
