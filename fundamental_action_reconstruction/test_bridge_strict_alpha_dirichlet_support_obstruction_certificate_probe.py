#!/usr/bin/env python3
"""Tests for the strict-alpha Dirichlet support obstruction certificate probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_dirichlet_support_obstruction_certificate_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_dirichlet_support_obstruction_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_dirichlet_support_obstruction_certificate_report.md"


class BridgeStrictAlphaDirichletSupportObstructionCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(
            self.payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_DIRICHLET_SUPPORT_OBSTRUCTION_CERTIFICATE_PROBE__NOT_A_THEOREM",
        )
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("This packet proves a support obstruction for that candidate; it does not derive a replacement support/phase-placement theorem.", self.payload["hard_limits"])
        self.assertIn("No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_target_and_previous_audit_replay(self) -> None:
        replay = self.payload["previous_dirichlet_audit_replay"]
        self.assertIn("STRICT_GATE_DIRICHLET_ACTION_AUDIT", replay["result_kind"])
        self.assertEqual(replay["global_best_ledger"], [2, 2, 2, 1, 1])
        self.assertTrue(replay["global_best_support_is_cyclic_contiguous"])

        target = self.payload["target_identity_replay"]
        self.assertEqual(target["q_power_product"], "256/243")
        self.assertLess(abs(target["eta_residual_vs_9_5"]), 1e-12)
        self.assertEqual(target["balanced_ledger"], [2, 2, 2, 1, 1])
        self.assertEqual(target["fifth_step_support_ordered"], [0, 7, 2, 9, 4])
        self.assertEqual(target["fifth_step_support_canonical"], [0, 2, 4, 7, 9])

    def test_indicator_support_scan_obstruction(self) -> None:
        indicator = self.payload["support_obstruction_scan"]["indicator_energy"]
        self.assertEqual(self.payload["support_obstruction_scan"]["support_count"], 792)
        self.assertEqual(indicator["minimizer_count"], 12)
        self.assertEqual(indicator["maximizer_count"], 12)
        self.assertTrue(indicator["contiguous_example_row"]["cyclic_contiguous"])
        self.assertFalse(indicator["fifth_step_row"]["cyclic_contiguous"])
        self.assertEqual(indicator["fifth_step_row"]["gap_signature"], [2, 2, 2, 3, 3])
        self.assertAlmostEqual(indicator["contiguous_example_row"]["energy"], 2.9296169079267154)
        self.assertAlmostEqual(indicator["fifth_step_row"]["energy"], 6.496452502177836)
        self.assertEqual(indicator["contiguous_example_row"]["energy"], indicator["min_energy"])
        self.assertEqual(indicator["fifth_step_row"]["energy"], indicator["max_energy"])

    def test_balanced_ledger_support_scan_obstruction(self) -> None:
        balanced = self.payload["support_obstruction_scan"]["balanced_ledger_min_assignment_energy"]
        self.assertEqual(self.payload["support_obstruction_scan"]["balanced_assignment_rows_scanned"], 7920)
        self.assertEqual(balanced["ledger"], [2, 2, 2, 1, 1])
        self.assertEqual(balanced["minimizer_count"], 12)
        self.assertEqual(balanced["maximizer_count"], 12)
        self.assertTrue(balanced["contiguous_example_row"]["cyclic_contiguous"])
        self.assertFalse(balanced["fifth_step_row"]["cyclic_contiguous"])
        self.assertAlmostEqual(balanced["contiguous_example_row"]["min_energy"], 8.066461680488555)
        self.assertAlmostEqual(balanced["fifth_step_row"]["min_energy"], 18.29375613429685)
        self.assertEqual(balanced["contiguous_example_row"]["min_energy"], balanced["min_energy"])
        self.assertEqual(balanced["fifth_step_row"]["min_energy"], balanced["max_of_support_minima"])

    def test_obstruction_summary_and_candidate_status(self) -> None:
        summary = self.payload["obstruction_summary"]
        self.assertTrue(summary["indicator_contiguous_minimum"])
        self.assertTrue(summary["indicator_fifth_step_is_maximum"])
        self.assertTrue(summary["balanced_contiguous_minimum"])
        self.assertTrue(summary["balanced_fifth_step_is_maximum_of_support_minima"])
        self.assertGreater(summary["energy_gap_indicator_fifth_minus_contiguous"], 0.0)
        self.assertGreater(summary["energy_gap_balanced_fifth_minus_contiguous"], 0.0)
        self.assertIn("clusters support", summary["verdict"])

        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("negative support certificate", candidate["content"])
        self.assertIn("all 792 five-node supports", candidate["why_this_is_more_proof_like"])
        self.assertIn("cannot be the missing strict-core support selector", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
