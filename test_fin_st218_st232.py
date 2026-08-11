#!/usr/bin/env python3
"""Lightweight live/serialized validation for FIN programs ST218--ST232."""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from fin_st01_st15_research import strict_operator


ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / "FIN_ST218_ST232_Results.json").read_text(encoding="utf-8"))
COUNT = 0


def check(condition: bool, message: str) -> None:
    global COUNT
    if not condition:
        raise AssertionError(message)
    COUNT += 1


def main() -> None:
    check(all(f"ST{k}" in DATA for k in range(218, 233)), "complete program range")
    check(DATA["metadata"]["programs"] == "ST218-ST232", "metadata range")

    for k in range(218, 233):
        row = DATA[f"ST{k}"]; path = ROOT / row["packet_file"]
        check(path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() == row["packet_sha256"], f"ST{k} packet hash")

    _, a, _ = strict_operator()
    check(DATA["ST218"]["strict_polynomial_algebra_dimension"] == 7 and DATA["ST218"]["full_commutant_complex_dimension"] == 22, "ST218 dimensions")
    check(max(x["commutator_Frobenius_norm"] for x in DATA["ST218"]["candidate_rows"]) < 1e-12, "ST218 candidates commute")
    check(np.linalg.norm(a @ expm(-0.37 * a) - expm(-0.37 * a) @ a) < 1e-12, "ST218 live functional calculus")

    check([x["minimal_atoms"] for x in DATA["ST219"]["rows"]] == [12, 7, 1, 1], "ST219 atom counts")
    check(DATA["ST219"]["rows"][0]["partition_lattice_size_Bell"] == 4213597 and DATA["ST219"]["rows"][1]["partition_lattice_size_Bell"] == 877, "ST219 Bell sizes")

    check(DATA["ST220"]["group_order"] == 24 and DATA["ST220"]["invariant_ordered_pair_orbits_and_covariance_dimension"] == 7, "ST220 group classification")
    check(DATA["ST220"]["fixed_vector_dimension_mean_zero_sector"] == 0 and DATA["ST220"]["mean_zero_isotropic_covariance_invariance_error"] == 0.0, "ST220 selector obstruction")

    st221 = DATA["ST221"]
    check(st221["certified_seed_count"] == 60 and st221["certified_steps"] == st221["attempted_steps"] == 360, "ST221 complete finite campaign")
    check(all(r["certificate"] and r["certificate"]["included"] for r in st221["rows"]), "ST221 all Krawczyk inclusions")
    check(st221["minimum_Jacobian_rank_margin"] > 0 and st221["minimum_tangent_alignment"] > 0.99, "ST221 rank and tangent margins")

    st222 = DATA["ST222"]
    check(sum(Fraction(x) for x in st222["exact_Pauli_frame_probabilities"]) == 1, "ST222 exact probabilities")
    check(st222["Bloch_determinant"] < 0 and abs(st222["primal_entanglement_fidelity"] - 0.41) < 1e-14, "ST222 negative orientation optimum")
    check(abs(st222["primal_dual_gap"]) < 1e-14 and st222["laboratory_offdiagonal_Frobenius_norm"] > 0.3, "ST222 certificate and non-diagonality")

    st223 = DATA["ST223"]
    check(st223["optimistic_interval_DP_lower_bound"] <= st223["constructive_schedule_cost_upper_bound"], "ST223 valid bracket ordering")
    check(st223["terminal_interval"][0] <= st223["constructive_final_probability"] <= st223["terminal_interval"][1], "ST223 feasible terminal")
    check(len(st223["constructive_schedule"]) == st223["steps"] and st223["bracket_width"] < 0.02, "ST223 schedule and resolution")

    check(max(r["invariant_transport_error"] for r in DATA["ST224"]["invertible_transport_rows"]) < 1e-12, "ST224 invertible transport")
    check(max(r["invariant_error"] for r in DATA["ST224"]["direction_rescaling_rows"]) < 1e-12, "ST224 projective direction")

    stage = DATA["ST225"]["stages"][0]
    check(stage["global_halfwidth"] == 0.0004 and stage["boxes"] == stage["passed"] == 28 ** 3, "ST225 full doubled cover")
    check(stage["minimum_margin"] > 9e-4 and DATA["ST225"]["largest_complete_halfwidth"] == 0.0004, "ST225 positive interval margin")

    risks = DATA["ST226"]["rows"]
    check(all(r["optimized_Le_Cam_risk_lower_bound"] >= 0 and r["bias_derivative_0_25_CR_lower_bound"] > 0 for r in risks), "ST226 nonnegative minimax bounds")
    check(all(risks[i + 1]["bias_derivative_0_25_CR_lower_bound"] < risks[i]["bias_derivative_0_25_CR_lower_bound"] for i in range(len(risks) - 1)), "ST226 sample monotonicity")

    st227 = DATA["ST227"]
    check(st227["minimum_transversality_gap"] > 0.4 and st227["constraint_error"] < 1e-12, "ST227 polar uniqueness premises")

    st228 = DATA["ST228"]
    check(st228["belief_cylinders"] == 65536 and abs(st228["width_reduction_factor"] - 0.75) < 1e-14, "ST228 exact belief cylinders")
    check(st228["new_bracket"]["asymptotic_rate_lower"] < st228["new_bracket"]["asymptotic_rate_upper"], "ST228 bracket ordering")

    st229 = DATA["ST229"]
    check(st229["final_global_SWAP_energy_commutator_norm"] == 0.0 and min(st229["individual_pulse_energy_commutator_norms"]) > 2, "ST229 net versus pointwise conservation")
    check(st229["random_state_reset_joint_error"] == 0.0 and "transferred" in st229["information_destination"], "ST229 reversible reset")

    st230 = DATA["ST230"]
    check(len(st230["rows"]) == 6 and st230["diffusion_over_resistance_range"][1] - st230["diffusion_over_resistance_range"][0] > 0.1, "ST230 nonproportional metrics")
    check(all(all(r[k] > 0 for k in ["diffusion_t1", "sqrt_effective_resistance", "Hellinger_Fisher_chord", "first_projector_chord"]) for r in st230["rows"]), "ST230 positive distances")

    st231 = DATA["ST231"]
    check(max(r["generator_intertwining_error"] for r in st231["rows"]) < 1e-13, "ST231 generator intertwining")
    check(max(r["distance_error"] for r in st231["rows"]) < 1e-13 and len({r["mu"] for r in st231["rows"]}) == 4, "ST231 exact geometry and nonuniqueness")

    check(DATA["ST232"]["status"] == "blocked_no_external_record" and not DATA["ST232"]["local_search_performed"], "ST232 evidence stop")
    check(len(DATA["recommended_next_programs"]) == 15 and DATA["recommended_next_programs"][0]["id"] == "ST233", "recommendations")
    boundary = DATA["epistemic_boundary"]
    check("QW-2191" in boundary and "ToE closure" in boundary and "legacy-to-strict" in boundary, "epistemic guardrail")
    check("single-operator" in DATA["central_verdict"] and "arbitrary fiber rate" in DATA["central_verdict"], "central verdict")

    print(f"{COUNT}/{COUNT} tests passed")


if __name__ == "__main__":
    main()
