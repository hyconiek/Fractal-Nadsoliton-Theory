#!/usr/bin/env python3
"""Tests for the strict-alpha d5 selector symmetry orbit audit probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d5_selector_symmetry_orbit_audit_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d5_selector_symmetry_orbit_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_d5_selector_symmetry_orbit_audit_report.md"


class BridgeStrictAlphaD5SelectorSymmetryOrbitAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(
            self.payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_D5_SELECTOR_SYMMETRY_ORBIT_AUDIT_PROBE__NOT_A_THEOREM",
        )
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("The source/orientation anchor is confirmed as extra selector data, not a strict-core consequence of D12 symmetry.", self.payload["hard_limits"])
        self.assertIn("No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged; this packet quantifies the remaining free orbit.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_previous_anchor_and_target_replay(self) -> None:
        replay = self.payload["previous_anchored_selector_replay"]
        self.assertIn("D5_ANCHORED_ENDPOINT_ASSIGNMENT_SELECTOR", replay["result_kind"])
        self.assertTrue(replay["all_forward_assignments_unique"])
        self.assertEqual(replay["forward_selected_assignment_set"], [[2, 2, 2, 1, 1]])

        target = self.payload["target_identity_replay"]
        self.assertEqual(target["q_power_product"], "256/243")
        self.assertLess(abs(target["eta_residual_vs_9_5"]), 1e-12)
        self.assertEqual(target["forward_assignment"], [2, 2, 2, 1, 1])

    def test_anchored_configuration_enumeration(self) -> None:
        enum = self.payload["anchored_configuration_enumeration"]
        self.assertEqual(enum["row_count"], 24)
        self.assertEqual(enum["unique_value_configuration_count"], 24)
        self.assertEqual(enum["unique_unoriented_support_count"], 12)
        self.assertEqual(enum["source_count"], 12)
        self.assertEqual(enum["orientation_count_per_source"], 2)
        self.assertEqual(len(enum["rows"]), 24)

    def test_dihedral_orbit_and_stabilizer(self) -> None:
        orbit = self.payload["dihedral_orbit_audit"]
        self.assertEqual(orbit["dihedral_group_order"], 24)
        self.assertEqual(orbit["value_configuration_orbit_size"], 24)
        self.assertEqual(orbit["value_configuration_stabilizer_size"], 1)
        self.assertEqual(orbit["value_configuration_stabilizer"], [{"shift": 0, "reflect": False}])
        self.assertEqual(orbit["support_orbit_size_without_orientation_assignment"], 12)
        self.assertEqual(len(orbit["support_orbit"]), 12)

    def test_symmetry_verdict_and_candidate_status(self) -> None:
        verdict = self.payload["symmetry_verdict"]
        self.assertTrue(verdict["anchored_rows_equal_full_dihedral_orbit"])
        self.assertTrue(verdict["value_configuration_has_trivial_stabilizer"])
        self.assertTrue(verdict["unoriented_support_orbit_has_12_members"])
        self.assertIn("free D12 orbit", verdict["interpretation"])

        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("free D12 orbit", candidate["content"])
        self.assertIn("D12 orbit size and stabilizer", candidate["why_this_is_more_proof_like"])
        self.assertIn("QW-2191 is not discharged", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
