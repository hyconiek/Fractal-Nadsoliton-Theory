#!/usr/bin/env python3
"""Tests for the scratch phase-normalized measure-transport bridge probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_phase_normalized_measure_transport_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_phase_normalized_measure_transport_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_phase_normalized_measure_transport_report.md"


class BridgePhaseNormalizedMeasureTransportProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(
            self.payload["result_kind"],
            "SCRATCH_PHASE_NORMALIZED_MEASURE_TRANSPORT_PROBE__NOT_A_THEOREM",
        )
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])
        self.assertFalse(self.payload["affine_failure"]["affine_budget_measure_transport_acceptable"])

    def test_exact_measure_balance_and_unit_budgets(self) -> None:
        balance = self.payload["measure_balance"]
        self.assertLess(balance["max_abs_balance_residual"], 1e-12)
        self.assertLess(abs(balance["strict_exact_residual_vs_1"]), 1e-12)
        self.assertLess(abs(balance["legacy_exact_residual_vs_1"]), 1e-12)
        self.assertLess(abs(balance["strict_trapezoid_residual_vs_1"]), 1e-4)
        self.assertLess(abs(balance["legacy_trapezoid_residual_vs_1"]), 1e-4)
        self.assertLess(abs(balance["transported_trapezoid_residual_vs_1"]), 1e-4)

    def test_affine_budget_is_discriminated(self) -> None:
        affine = self.payload["affine_failure"]
        self.assertGreater(affine["affine_balance_max_abs_residual"], 1.0)
        self.assertGreater(affine["affine_balance_mean_abs_residual"], 0.1)

    def test_candidate_and_replays(self) -> None:
        candidate = self.payload["candidate_ontological_reading"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("transport", candidate["content"])
        replays = self.payload["upstream_replay"]
        self.assertTrue(replays["phase_budget_candidate_supported"])
        self.assertTrue(replays["horizon_candidate_supported"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
