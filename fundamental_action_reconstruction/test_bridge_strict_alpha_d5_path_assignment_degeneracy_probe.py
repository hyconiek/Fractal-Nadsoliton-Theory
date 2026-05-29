#!/usr/bin/env python3
"""Tests for the strict-alpha d5 path assignment degeneracy probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d5_path_assignment_degeneracy_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d5_path_assignment_degeneracy_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d5_path_assignment_degeneracy_report.md"


class BridgeStrictAlphaD5PathAssignmentDegeneracyProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(
            self.payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_D5_PATH_ASSIGNMENT_DEGENERACY_PROBE__NOT_A_THEOREM",
        )
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("The path-smooth assignment selector is a conditional model-internal premise, not a derived strict nadsoliton action theorem.", self.payload["hard_limits"])
        self.assertIn("Distance-5 support plus path-smooth assignment still leaves endpoint/translate/orientation degeneracy.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_bandpass_and_target_replay(self) -> None:
        replay = self.payload["previous_bandpass_replay"]
        self.assertIn("BANDPASS_D5_SUPPORT_SELECTOR", replay["result_kind"])
        self.assertTrue(replay["maximizers_equal_fifth_step_orbit"])
        self.assertEqual(replay["max_h5"], 4)
        self.assertEqual(replay["maximizer_count"], 12)

        target = self.payload["target_identity_replay"]
        self.assertEqual(target["q_power_product"], "256/243")
        self.assertLess(abs(target["eta_residual_vs_9_5"]), 1e-12)
        self.assertEqual(target["balanced_ledger"], [2, 2, 2, 1, 1])

    def test_path_assignment_model_for_reference_support(self) -> None:
        model = self.payload["path_assignment_model"]
        self.assertEqual(model["ordered_support_example"], [0, 7, 2, 9, 4])
        self.assertTrue(model["all_adjacent_distances_are_5"])
        self.assertEqual(model["energy"], "E_path(a)=sum_{r=0..3}(a_{r+1}-a_r)^2")
        self.assertTrue(model["energy_equals_transition_count_for_values_1_and_2"])
        self.assertEqual(model["transition_histogram_for_one_ordered_support"], {"1": 2, "2": 3, "3": 4, "4": 1})
        self.assertEqual(len(model["reference_rows"]), 10)

    def test_orbit_assignment_scan_degeneracy(self) -> None:
        scan = self.payload["orbit_assignment_scan"]
        self.assertEqual(scan["ordered_support_count"], 24)
        self.assertEqual(scan["assignments_per_ordered_support"], 10)
        self.assertEqual(scan["total_assignment_rows"], 240)
        self.assertEqual(scan["min_path_variation_energy"], 1)
        self.assertEqual(scan["max_path_variation_energy"], 4)
        self.assertEqual(scan["minimizer_count"], 48)
        self.assertEqual(scan["maximizer_count"], 24)
        self.assertEqual(scan["unique_minimizer_high_node_sets_count"], 12)
        self.assertEqual(len(scan["unique_minimizer_high_node_sets"]), 12)

    def test_reference_minimizers_and_maximizer(self) -> None:
        scan = self.payload["orbit_assignment_scan"]
        self.assertEqual(
            [row["assignment"] for row in scan["reference_minimizers"]],
            [[1, 1, 2, 2, 2], [2, 2, 2, 1, 1]],
        )
        self.assertEqual(
            [row["high_value_nodes"] for row in scan["reference_minimizers"]],
            [[2, 4, 9], [0, 2, 7]],
        )
        self.assertEqual([row["assignment"] for row in scan["reference_maximizers"]], [[2, 1, 2, 1, 2]])
        self.assertEqual(scan["reference_maximizers"][0]["path_variation_energy"], 4)

    def test_conditional_schema_and_candidate_status(self) -> None:
        schema = self.payload["conditional_assignment_schema"]
        self.assertIn("distance-5 band-pass support selector", schema["support_premise"])
        self.assertIn("(2,2,2,1,1)", schema["ledger_premise"])
        self.assertIn("minimize E_path", schema["path_smooth_assignment_premise"])
        self.assertIn("22211 or 11222", schema["conditional_output"])
        self.assertIn("endpoint side", schema["remaining_obstruction"])

        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("controlled degeneracy", candidate["content"])
        self.assertIn("240-row orbit scan", candidate["why_this_is_more_proof_like"])
        self.assertIn("No strict theorem derives", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
