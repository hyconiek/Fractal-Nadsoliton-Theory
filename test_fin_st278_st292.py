#!/usr/bin/env python3
"""Independent live checks for FIN ST278--ST292."""

from __future__ import annotations

import hashlib
import json
import math
import unittest
from pathlib import Path

import numpy as np
import sympy as sp

from fin_st01_st15_research import strict_operator
from fin_st278_st292_research import pauli
from validate_fin_st276_record import validate


ROOT = Path(__file__).resolve().parent
D = json.loads((ROOT / "FIN_ST278_ST292_Results.json").read_text())


class TestST278ST292(unittest.TestCase):
    def test_packets_and_hashes(self):
        for k in range(278, 293):
            row = D[f"ST{k}"]
            path = ROOT / row["packet_file"]
            self.assertTrue(path.exists())
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), row["packet_sha256"])

    def test_st278_exact_coefficient_and_threshold(self):
        self.assertEqual(sp.Rational(D["ST278"]["exact_degree_12_coefficient"]), sp.Rational(7776, 11))
        _, a, _ = strict_operator()
        lam = np.linalg.eigvalsh(a)[1]
        self.assertAlmostEqual(D["ST278"]["critical_coupling_g"], 12 / lam, places=11)
        self.assertGreater(D["ST278"]["unit_coupling_curvature"], 0)

    def test_st279_generators_enumerate_covariants(self):
        for row in D["ST279"]["enumeration"]:
            self.assertEqual((row["p"] - row["q"] - 1) % 12, 0)
        self.assertEqual(D["ST279"]["Hilbert_basis"], ["z", "conjugate(z)^11"])

    def test_st280_positive_rank_boxes(self):
        bounds = D["ST280"]["certified_interval_rank_lower_bounds"]
        self.assertEqual(D["ST280"]["certified_three_box_steps"], [21, 22, 23])
        self.assertTrue(all(x > 0 for x in bounds))
        self.assertLess(bounds[1], bounds[0])
        self.assertLess(bounds[1], bounds[2])

    def test_st281_formula_examples(self):
        for row in D["ST281"]["examples"]:
            self.assertAlmostEqual(row["optimal_entanglement_fidelity"],
                                   (1 + row["trace_norm"] / 2) / 4, places=14)

    def test_st282_positive_rationalized_margin(self):
        self.assertGreater(D["ST282"]["minimum_Krawczyk_margin"], 0)
        self.assertLess(D["ST282"]["preconditioner_inverse_defect_inf_norm"], 1e-9)

    def test_st283_exact_circle_counterexample(self):
        self.assertAlmostEqual(D["ST283"]["second_fundamental_form_norm"], 1.0)
        self.assertLess(D["ST283"]["tangential_target_derivative"], 1e-14)

    def test_st284_complete_cover(self):
        self.assertTrue(D["ST284"]["complete_cover"])
        self.assertEqual(D["ST284"]["final_failed_boxes"], 0)
        self.assertEqual(D["ST284"]["adaptive_levels"][0]["tested_boxes"], 64000)
        self.assertGreater(D["ST284"]["adaptive_levels"][-1]["minimum_margin"], 0)

    def test_st286_live_clifford_lift(self):
        sx, sy, sz = pauli()
        phi, a, b = 0.371, 3 / 5, 4 / 5
        u = math.cos(phi) * sx + math.sin(phi) * sy
        X, Z = a * sz + b * u, b * sz - a * u
        self.assertLess(np.linalg.norm(X @ X - np.eye(2)), 1e-14)
        self.assertLess(np.linalg.norm(Z @ Z - np.eye(2)), 1e-14)
        self.assertLess(np.linalg.norm(X @ Z + Z @ X), 1e-14)

    def test_st287_endpoint_derivative_signs(self):
        left, right = D["ST287"]["derivative_audits"]
        self.assertLess(left["derivative_global_upper"], 0)
        self.assertGreater(right["derivative_global_lower"], 0)

    def test_st288_controller(self):
        self.assertLess(D["ST288"]["unit_time_error_from_controlled_SWAP"], 1e-12)
        self.assertLess(D["ST288"]["controller_commutator_norm"], 1e-14)

    def test_st289_modulus_not_selected(self):
        self.assertLess(D["ST289"]["associator_defect"], 1e-14)
        self.assertGreater(D["ST289"]["unequal_rate_fiber_swap_defect"], 0.1)
        self.assertLess(D["ST289"]["equal_rate_fiber_swap_defect"], 1e-14)

    def test_st290_tradeoff_monotonicity(self):
        trace = np.array(D["ST290"]["trace_budget_at_t1"])
        end = np.array(D["ST290"]["guaranteed_plateau_endpoint"])
        self.assertTrue(np.all(np.diff(trace) < 0))
        self.assertTrue(np.all(np.diff(end) < 0))
        self.assertTrue(all(x["numeric_heat_trace_t1"] < math.inf for x in D["ST290"]["soft_power_audit"]))

    def test_st291_validator_fails_closed(self):
        record = {"schema_version": "FIN-ST276-1", "provider": "same", "registrar": "same",
                  "analyst": "same", "holdout_frozen": False, "calibrated_time": -1,
                  "run_id": "x", "calibration_hash": "bad", "protocol_hash": "bad",
                  "events": [], "synthetic": True}
        errors = validate(record, True, True)
        self.assertGreaterEqual(len(errors), 5)
        self.assertTrue(any("synthetic" in x for x in errors))

    def test_st292_blocked(self):
        self.assertEqual(D["ST292"]["independent_record_count"], 0)
        self.assertTrue(D["ST292"]["status"].startswith("blocked"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
