#!/usr/bin/env python3
"""Tests for the strict-alpha phase/Fourier selector discriminator probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_phase_fourier_selector_discriminator_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_phase_fourier_selector_discriminator_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_phase_fourier_selector_discriminator_report.md"


class BridgeStrictAlphaPhaseFourierSelectorDiscriminatorProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(
            self.payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_PHASE_FOURIER_SELECTOR_DISCRIMINATOR_PROBE__NOT_A_THEOREM",
        )
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No theorem derives a strict-core phase/Fourier selector from nadsoliton geometry.", self.payload["hard_limits"])
        self.assertIn("Naive highest total Fourier resonance selects (4,1,1,1,1), not (2,2,2,1,1).", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_target_identity_replay(self) -> None:
        target = self.payload["target_identity_replay"]
        self.assertEqual(target["q_radical_power_5"], "256/243")
        self.assertEqual(target["all_canonical_ledgers_share_product"], "2^8/3^5 = 256/243")
        self.assertAlmostEqual(target["eta_target"], 9.0 / 5.0)
        self.assertIn("PHASE_RESONANCE_LIMMATIC_AUDIT", target["phase_audit_result_kind"])
        self.assertIn("ENTROPY_SELECTOR_DISCRIMINATOR", target["entropy_audit_result_kind"])

    def test_parseval_discriminator_selects_minimum_not_highest(self) -> None:
        parseval = self.payload["parseval_discriminator"]
        self.assertEqual(parseval["identity"], "total_non_dc_power = N*sum(e_i^2) - (sum e_i)^2")
        self.assertEqual(parseval["minimum_fourier_ripple_winner"]["ledger"], [2, 2, 2, 1, 1])
        self.assertEqual(parseval["maximum_fourier_ripple_winner"]["ledger"], [4, 1, 1, 1, 1])
        self.assertEqual(parseval["z12_zero_padded_minimum_fourier_ripple_winner"]["ledger"], [2, 2, 2, 1, 1])
        self.assertEqual(parseval["z12_zero_padded_maximum_fourier_ripple_winner"]["ledger"], [4, 1, 1, 1, 1])
        self.assertIn("MINIMUM_RIPPLE_NOT_HIGHEST", parseval["verdict"])

        rows = {tuple(row["ledger"]): row for row in parseval["ledger_rows"]}
        self.assertEqual(rows[(4, 1, 1, 1, 1)]["sum_squares"], 20)
        self.assertEqual(rows[(3, 2, 1, 1, 1)]["sum_squares"], 16)
        self.assertEqual(rows[(2, 2, 2, 1, 1)]["sum_squares"], 14)
        self.assertEqual(rows[(4, 1, 1, 1, 1)]["five_branch_total_non_dc_power_parseval"], 36)
        self.assertEqual(rows[(3, 2, 1, 1, 1)]["five_branch_total_non_dc_power_parseval"], 16)
        self.assertEqual(rows[(2, 2, 2, 1, 1)]["five_branch_total_non_dc_power_parseval"], 6)
        for row in rows.values():
            self.assertEqual(row["product_fraction"], "256/243")
            self.assertLess(abs(row["eta_residual_vs_9_5"]), 1e-12)

    def test_single_mode_scans_are_placement_sensitive(self) -> None:
        fixed = self.payload["fixed_fifth_window_single_mode_scan"]
        self.assertIn("phase-labelled branch slots", fixed["interpretation"])
        self.assertIn("SINGLE_MODE_CRITERIA_DEPEND_ON_ASSIGNMENT", fixed["verdict"])
        self.assertEqual([row["fifth_window_nodes"] for row in fixed["rows"]], [[0, 7, 2, 9, 4]] * 3)

        bounded = self.payload["bounded_z12_placement_scan"]
        self.assertTrue(bounded["all_ledgers_can_reach_dc_sized_single_mode_power_64"])
        self.assertIn("PLACEMENT_DOMINATED", bounded["verdict"])
        rows = {tuple(row["ledger"]): row for row in bounded["rows"]}
        self.assertEqual(rows[(4, 1, 1, 1, 1)]["rows_scanned"], 3960)
        self.assertEqual(rows[(3, 2, 1, 1, 1)]["rows_scanned"], 15840)
        self.assertEqual(rows[(2, 2, 2, 1, 1)]["rows_scanned"], 7920)
        for row in rows.values():
            self.assertAlmostEqual(row["max_possible_max_single_non_dc_power"], 64.0)
            self.assertGreater(row["distinct_rounded_max_single_values"], 1)

    def test_candidate_supported_but_not_closure(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("minimum ripple", candidate["content"])
        self.assertIn("Parseval identity", candidate["why_this_is_more_proof_like"])
        self.assertIn("No strict theorem derives", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
