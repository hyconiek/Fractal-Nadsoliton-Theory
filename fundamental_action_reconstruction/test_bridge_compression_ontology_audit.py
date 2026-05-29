#!/usr/bin/env python3
"""Tests for the scratch compression ontology audit."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_compression_ontology_audit.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_compression_ontology_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_compression_ontology_audit_report.md"


class BridgeCompressionOntologyAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_no_false_pass(self) -> None:
        self.assertEqual(self.payload["result_kind"], "SCRATCH_COMPRESSION_ONTOLOGY_AUDIT__NOT_A_THEOREM")
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No legacy physical-role transfer is licensed.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_strict_has_stronger_operational_compression(self) -> None:
        compression = self.payload["compression_characteristic"]
        legacy_d11 = compression["missing_in_legacy_reading"]["legacy_denominator_at_z12_max_d11"]
        strict_d11 = compression["present_in_strict_operational_form"]["strict_denominator_at_z12_max_d11"]
        ratio = compression["present_in_strict_operational_form"]["strict_vs_legacy_denominator_ratio_at_d11"]
        self.assertGreater(legacy_d11, 0.8)
        self.assertLess(strict_d11, 0.02)
        self.assertLess(ratio, 0.02)

    def test_measure_transport_replay_and_density_moments(self) -> None:
        replay = self.payload["upstream_replay"]
        moments = self.payload["strict_loss_density_moments"]
        self.assertTrue(replay["measure_transport_candidate_supported"])
        self.assertLess(replay["measure_balance_max_abs_residual"], 1e-12)
        self.assertTrue(replay["phase_budget_candidate_supported"])
        self.assertLess(abs(moments["trapezoid_norm"] - 1.0), 1e-4)
        self.assertGreater(moments["left_half_mass"], 0.9)
        self.assertLess(moments["mean_v"], 0.25)

    def test_final_form_assessment_is_honest(self) -> None:
        assessment = self.payload["strict_final_form_assessment"]
        self.assertTrue(assessment["is_strict_plausible_terminal_shape_candidate"])
        self.assertFalse(assessment["is_strict_proven_final_nadsoliton_formula"])
        self.assertGreaterEqual(len(assessment["blocking_obligations"]), 4)
        answers = " ".join(self.payload["honest_answer_to_user_question"])
        self.assertIn("missing explicit characteristic", answers)
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
