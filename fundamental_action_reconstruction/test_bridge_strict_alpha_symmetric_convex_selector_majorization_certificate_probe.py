#!/usr/bin/env python3
"""Tests for the strict-alpha symmetric-convex majorization certificate probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_symmetric_convex_selector_majorization_certificate_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_symmetric_convex_selector_majorization_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_symmetric_convex_selector_majorization_certificate_report.md"


class BridgeStrictAlphaSymmetricConvexSelectorMajorizationCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(
            self.payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_SYMMETRIC_CONVEX_SELECTOR_MAJORIZATION_CERTIFICATE_PROBE__NOT_A_THEOREM",
        )
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No theorem derives a symmetric strictly convex branch action from strict nadsoliton geometry.", self.payload["hard_limits"])
        self.assertIn("This certificate is conditional on a Schur-convex selector/action premise; it is not a strict-core selector source.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_conditional_theorem_statement(self) -> None:
        theorem = self.payload["conditional_theorem_statement"]
        self.assertIn("symmetric strictly convex separable branch action, or any Schur-convex branch action", theorem["premises"])
        self.assertIn("majorizes the floor/ceil balanced ledger", theorem["majorization_claim"])
        self.assertIn("strictly decreases every checked convex power energy", theorem["smoothing_claim"])
        self.assertIn("symmetric strictly convex selector", theorem["conclusion"])
        self.assertIn("strict nadsoliton action source is not derived", theorem["physical_status"])

    def test_target_majorization_certificate(self) -> None:
        target = self.payload["target_eta_9_5_certificate"]
        self.assertEqual(target["q_power_product"], "256/243")
        self.assertLess(abs(target["eta_residual_vs_9_5"]), 1e-12)
        self.assertEqual(target["branch_count"], 5)
        self.assertEqual(target["total_binary_exponent"], 8)
        self.assertEqual(target["balanced_ledger"], [2, 2, 2, 1, 1])
        self.assertTrue(target["all_ledgers_majorize_balanced"])
        self.assertTrue(target["balanced_majorizes_only_itself"])
        self.assertTrue(target["all_checked_powers_select_balanced"])
        for power in [2, 3, 4, 5, 6]:
            self.assertEqual(target["unique_minimizers_by_power"][f"p_{power}"]["minimizers"], [[2, 2, 2, 1, 1]])

    def test_power_energies_and_smoothing_paths(self) -> None:
        rows = {tuple(row["ledger"]): row for row in self.payload["target_eta_9_5_certificate"]["rows"]}
        self.assertEqual(rows[(2, 2, 2, 1, 1)]["power_energies"]["p_2"], 14)
        self.assertEqual(rows[(3, 2, 1, 1, 1)]["power_energies"]["p_2"], 16)
        self.assertEqual(rows[(4, 1, 1, 1, 1)]["power_energies"]["p_2"], 20)
        self.assertEqual(rows[(2, 2, 2, 1, 1)]["power_energies"]["p_3"], 26)
        self.assertEqual(rows[(3, 2, 1, 1, 1)]["power_energies"]["p_3"], 38)
        self.assertEqual(rows[(4, 1, 1, 1, 1)]["power_energies"]["p_3"], 68)
        for row in rows.values():
            self.assertEqual(row["product_fraction"], "256/243")
            self.assertLess(abs(row["eta_residual_vs_9_5"]), 1e-12)
            self.assertEqual(row["smoothing_terminal"], [2, 2, 2, 1, 1])

        paths = self.payload["target_eta_9_5_certificate"]["sample_smoothing_paths"]
        self.assertEqual(paths["from_4_1_1_1_1"][0]["energy_drops"], {"p_2": 4, "p_3": 30, "p_4": 160, "p_5": 750, "p_6": 3304})
        self.assertEqual(paths["from_3_2_1_1_1"][0]["energy_drops"], {"p_2": 2, "p_3": 12, "p_4": 50, "p_5": 180, "p_6": 602})
        self.assertEqual(paths["from_2_2_2_1_1"], [])

    def test_bounded_majorization_scan(self) -> None:
        scan = self.payload["bounded_majorization_scan"]
        self.assertEqual(scan["scan_limits"], {"max_branch_count": 9, "max_total_exponent": 24, "powers_checked": [2, 3, 4, 5, 6]})
        self.assertEqual(scan["cases_checked"], 180)
        self.assertEqual(scan["canonical_ledgers_checked"], 5583)
        self.assertEqual(scan["max_partition_count"], 201)
        self.assertEqual(scan["max_partition_case"], {"branch_count": 7, "total": 24})
        self.assertTrue(scan["all_cases_passed"])
        self.assertEqual(scan["failures"], [])
        self.assertGreaterEqual(len(scan["sample_rows"]), 6)

    def test_candidate_supported_but_not_closure(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("symmetric strictly convex", candidate["content"])
        self.assertIn("majorization theorem class", candidate["why_this_is_more_proof_like"])
        self.assertIn("does not derive", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
