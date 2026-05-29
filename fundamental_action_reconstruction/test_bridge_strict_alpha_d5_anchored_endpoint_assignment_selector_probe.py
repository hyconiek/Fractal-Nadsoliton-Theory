#!/usr/bin/env python3
"""Tests for the strict-alpha d5 anchored endpoint assignment selector probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d5_anchored_endpoint_assignment_selector_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d5_anchored_endpoint_assignment_selector_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d5_anchored_endpoint_assignment_selector_report.md"


class BridgeStrictAlphaD5AnchoredEndpointAssignmentSelectorProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(
            self.payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_D5_ANCHORED_ENDPOINT_ASSIGNMENT_SELECTOR_PROBE__NOT_A_THEOREM",
        )
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("The anchored endpoint assignment selector is conditional on explicit source, orientation, and endpoint-bias premises.", self.payload["hard_limits"])
        self.assertIn("No theorem derives the endpoint-bias/source-moment action from strict nadsoliton geometry.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_previous_degeneracy_and_target_replay(self) -> None:
        replay = self.payload["previous_path_degeneracy_replay"]
        self.assertIn("D5_PATH_ASSIGNMENT_DEGENERACY", replay["result_kind"])
        self.assertEqual(replay["minimizer_count"], 48)
        self.assertEqual(replay["unique_minimizer_high_node_sets_count"], 12)

        target = self.payload["target_identity_replay"]
        self.assertEqual(target["q_power_product"], "256/243")
        self.assertLess(abs(target["eta_residual_vs_9_5"]), 1e-12)
        self.assertEqual(target["balanced_ledger"], [2, 2, 2, 1, 1])

    def test_selector_definition_and_reference_certificate(self) -> None:
        definition = self.payload["anchored_selector_definition"]
        self.assertIn("source + r*orientation*5", definition["support_rule"])
        self.assertEqual(definition["path_energy"], "E_path(a)=sum_r(a_{r+1}-a_r)^2")
        self.assertEqual(definition["source_moment"], "M_source(a)=sum_r r*a_r")
        self.assertIn("lexicographically minimize", definition["forward_endpoint_selector"])

        ref = self.payload["reference_source_orientation_certificate"]
        self.assertEqual(ref["source"], 0)
        self.assertEqual(ref["orientation"], -1)
        self.assertEqual(ref["ordered_support"], [0, 7, 2, 9, 4])
        self.assertTrue(ref["all_adjacent_distances_are_5"])
        self.assertEqual(ref["row_count"], 10)
        self.assertEqual(ref["path_min_energy"], 1)
        self.assertEqual(ref["path_minimizer_count_before_endpoint_bias"], 2)
        self.assertTrue(ref["forward_assignment_unique"])
        self.assertTrue(ref["reverse_assignment_unique"])
        self.assertEqual(ref["forward_endpoint_selected_row"]["assignment"], [2, 2, 2, 1, 1])
        self.assertEqual(ref["reverse_endpoint_selected_row"]["assignment"], [1, 1, 2, 2, 2])
        self.assertEqual(ref["forward_endpoint_selected_row"]["lex_forward_key"], [1, 13])
        self.assertEqual(ref["reverse_endpoint_selected_row"]["lex_reverse_key"], [1, 13])

    def test_orbit_scan_closes_assignment_given_anchor(self) -> None:
        scan = self.payload["orbit_scan"]
        self.assertEqual(scan["source_count"], 12)
        self.assertEqual(scan["orientation_count_per_source"], 2)
        self.assertEqual(scan["ordered_support_count"], 24)
        self.assertEqual(scan["assignments_per_ordered_support"], 10)
        self.assertTrue(scan["all_forward_assignments_unique"])
        self.assertTrue(scan["all_reverse_assignments_unique"])
        self.assertEqual(scan["forward_selected_assignment_set"], [[2, 2, 2, 1, 1]])
        self.assertEqual(scan["reverse_selected_assignment_set"], [[1, 1, 2, 2, 2]])
        self.assertEqual(len(scan["rows_summary"]), 24)

    def test_conditional_schema_and_candidate_status(self) -> None:
        schema = self.payload["conditional_selector_schema"]
        self.assertIn("explicit source and orientation", schema["support_premise"])
        self.assertIn("endpoint bias", schema["endpoint_bias_premise"])
        self.assertIn("22211", schema["conditional_output"])
        self.assertIn("strict origin", schema["remaining_obstruction"])

        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("remove the final finite assignment degeneracy", candidate["content"])
        self.assertIn("source moment", candidate["why_this_is_more_proof_like"])
        self.assertIn("remain premises", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
