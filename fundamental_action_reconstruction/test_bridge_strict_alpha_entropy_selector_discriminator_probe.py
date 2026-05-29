#!/usr/bin/env python3
"""Tests for the strict-alpha entropy selector discriminator scratch probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_entropy_selector_discriminator_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_entropy_selector_discriminator_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_entropy_selector_discriminator_report.md"


class BridgeStrictAlphaEntropySelectorDiscriminatorProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(self.payload["result_kind"], "SCRATCH_STRICT_ALPHA_ENTROPY_SELECTOR_DISCRIMINATOR_PROBE__NOT_A_THEOREM")
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No legacy physical-role transfer is licensed.", self.payload["hard_limits"])
        self.assertIn("No theorem derives fixed-labelled branch entropy as the strict-core selector.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_candidate_space_rows_share_exact_product(self) -> None:
        self.assertEqual(self.payload["candidate_space"]["canonical_positive_ledgers"], [[4, 1, 1, 1, 1], [3, 2, 1, 1, 1], [2, 2, 2, 1, 1]])
        for row in self.payload["selector_rows"]:
            self.assertEqual(row["product_fraction"], "256/243")
            self.assertLess(abs(row["eta_residual_vs_9_5"]), 1e-12)

    def test_fixed_labelled_entropy_selects_balanced_ledger(self) -> None:
        selector = self.payload["fixed_labelled_branch_entropy_selector"]
        self.assertTrue(selector["entropy_winner_unique"])
        self.assertTrue(selector["labelled_multinomial_winner_unique"])
        self.assertEqual(selector["entropy_winner"]["ledger"], [2, 2, 2, 1, 1])
        self.assertEqual(selector["labelled_multinomial_winner"]["ledger"], [2, 2, 2, 1, 1])
        self.assertEqual(selector["labelled_multinomial_winner"]["labelled_multinomial_count"], 5040)
        self.assertTrue(selector["matches_balanced_ledger"])

    def test_unlabelled_orbit_aggregate_selects_different_ledger(self) -> None:
        discriminator = self.payload["unlabelled_orbit_aggregate_discriminator"]
        self.assertTrue(discriminator["orbit_aggregate_winner_unique"])
        self.assertEqual(discriminator["orbit_aggregate_winner"]["ledger"], [3, 2, 1, 1, 1])
        self.assertEqual(discriminator["orbit_aggregate_winner"]["orbit_aggregate_count"], 67200)
        self.assertFalse(discriminator["matches_balanced_ledger"])
        self.assertIn("not the balanced ledger", discriminator["obstruction"])

    def test_candidate_supported_but_selector_material(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("fixed-labelled Shannon/multinomial selector", candidate["content"])
        self.assertIn("selector convention", candidate["why_this_is_more_proof_like"])
        self.assertIn("No strict theorem", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
