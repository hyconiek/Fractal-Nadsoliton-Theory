#!/usr/bin/env python3
"""Tests for the strict-alpha support distance-histogram polarity certificate probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_support_distance_histogram_polarity_certificate_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_support_distance_histogram_polarity_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_support_distance_histogram_polarity_certificate_report.md"


class BridgeStrictAlphaSupportDistanceHistogramPolarityCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(
            self.payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_SUPPORT_DISTANCE_HISTOGRAM_POLARITY_CERTIFICATE_PROBE__NOT_A_THEOREM",
        )
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("This packet explains a support obstruction for the strict-gate Dirichlet candidate; it does not derive a replacement support/phase-placement theorem.", self.payload["hard_limits"])
        self.assertIn("The polarity condition is a model-internal support-functional statement, not a strict nadsoliton action theorem.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_target_and_previous_obstruction_replay(self) -> None:
        replay = self.payload["previous_support_obstruction_replay"]
        self.assertIn("DIRICHLET_SUPPORT_OBSTRUCTION", replay["result_kind"])
        self.assertTrue(replay["indicator_fifth_step_is_maximum"])
        self.assertTrue(replay["balanced_fifth_step_is_maximum_of_support_minima"])

        target = self.payload["target_identity_replay"]
        self.assertEqual(target["q_power_product"], "256/243")
        self.assertLess(abs(target["eta_residual_vs_9_5"]), 1e-12)

    def test_linear_histogram_identity(self) -> None:
        identity = self.payload["linear_histogram_identity"]
        self.assertEqual(identity["support_size"], 5)
        self.assertEqual(identity["support_count"], 792)
        self.assertEqual(identity["histogram_class_count"], 35)
        self.assertEqual(identity["indicator_energy_formula"], "E_ind(S)=5*T-2*sum_d h_d(S)*w_d")
        weights = identity["distance_weights"]
        self.assertGreater(weights["w_1"], weights["w_5"])
        self.assertGreater(identity["total_node_weight_T"], 0.0)

    def test_contiguous_vs_fifth_exact_gap(self) -> None:
        cert = self.payload["contiguous_vs_fifth_certificate"]
        self.assertEqual(cert["contiguous_support"], [0, 1, 2, 3, 4])
        self.assertEqual(cert["contiguous_histogram"], [4, 3, 2, 1, 0, 0])
        self.assertTrue(cert["contiguous_is_global_indicator_minimum"])
        self.assertEqual(cert["fifth_step_support_ordered"], [0, 7, 2, 9, 4])
        self.assertEqual(cert["fifth_step_support_canonical"], [0, 2, 4, 7, 9])
        self.assertEqual(cert["fifth_step_histogram"], [0, 3, 2, 1, 4, 0])
        self.assertTrue(cert["fifth_step_is_global_indicator_maximum"])
        self.assertEqual(cert["histogram_difference_fifth_minus_contiguous"], [-4, 0, 0, 0, 4, 0])
        self.assertEqual(cert["exact_two_orbit_gap_formula"], "E_ind(fifth)-E_ind(contiguous)=8*(w_1-w_5)")
        self.assertAlmostEqual(cert["computed_gap"], cert["formula_gap"])
        self.assertAlmostEqual(cert["computed_gap"], 3.5668355942511214)
        self.assertIn("requires w_5>w_1", cert["polarity_condition_for_plain_minimization_to_prefer_fifth_over_contiguous"])

    def test_extremal_histogram_classes(self) -> None:
        extremal = self.payload["extremal_histogram_classes"]
        self.assertEqual(extremal["contiguous_histogram_support_count"], 12)
        self.assertEqual(extremal["fifth_step_histogram_support_count"], 12)
        self.assertEqual(len(extremal["minimizer_classes"]), 1)
        self.assertEqual(len(extremal["maximizer_classes"]), 1)
        self.assertEqual(extremal["minimizer_classes"][0]["distance_histogram"], [4, 3, 2, 1, 0, 0])
        self.assertEqual(extremal["maximizer_classes"][0]["distance_histogram"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(len(extremal["fifth_step_orbit_supports"]), 12)
        self.assertIn([0, 2, 4, 7, 9], extremal["fifth_step_orbit_supports"])

    def test_candidate_supported_but_not_closure(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("distance-histogram polarity problem", candidate["content"])
        self.assertIn("linear histogram functional", candidate["why_this_is_more_proof_like"])
        self.assertIn("does not derive a new support selector", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
