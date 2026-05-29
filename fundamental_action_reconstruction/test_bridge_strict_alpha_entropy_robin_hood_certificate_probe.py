#!/usr/bin/env python3
"""Tests for the strict-alpha entropy Robin-Hood certificate scratch probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_entropy_robin_hood_certificate_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_entropy_robin_hood_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_entropy_robin_hood_certificate_report.md"


class BridgeStrictAlphaEntropyRobinHoodCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(self.payload["result_kind"], "SCRATCH_STRICT_ALPHA_ENTROPY_ROBIN_HOOD_CERTIFICATE_PROBE__NOT_A_THEOREM")
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("No legacy physical-role transfer is licensed.", self.payload["hard_limits"])
        self.assertIn("No theorem derives fixed-labelled branch entropy as the strict-core selector.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_weight_law_is_conditional_not_selector_source(self) -> None:
        law = self.payload["fixed_labelled_weight_law"]
        self.assertEqual(law["weight"], "W(e_1,...,e_5)=8!/prod_i(e_i!)")
        self.assertIn("a/(b+1)", law["robin_hood_ratio"])
        self.assertTrue(law["derived_here_as_arithmetic_certificate"])
        self.assertFalse(law["strict_selector_source_derived_here"])

    def test_certificate_paths_increase_weight_and_entropy(self) -> None:
        summary = self.payload["certificate_summary"]
        self.assertEqual(summary["balanced_terminal"], [2, 2, 2, 1, 1])
        self.assertTrue(summary["source_entropy_winner_matches_terminal"])
        self.assertTrue(summary["all_nonterminal_ledgers_have_path"])
        self.assertTrue(summary["all_steps_increase_labelled_weight"])
        self.assertTrue(summary["all_steps_increase_entropy"])
        self.assertEqual(summary["terminal_labelled_weight"], 5040)

    def test_first_path_has_exact_ratios_and_weights(self) -> None:
        paths = self.payload["certificate_paths_to_balanced_ledger"]
        path = paths["(4, 1, 1, 1, 1)"]
        self.assertEqual([step["before"] for step in path], [[4, 1, 1, 1, 1], [3, 2, 1, 1, 1]])
        self.assertEqual([step["after"] for step in path], [[3, 2, 1, 1, 1], [2, 2, 2, 1, 1]])
        self.assertEqual([step["exact_labelled_weight_ratio_after_over_before"] for step in path], ["2/1", "3/2"])
        self.assertEqual([step["labelled_weight_after"] for step in path], [3360, 5040])
        self.assertTrue(all(step["entropy_increases"] for step in path))

    def test_candidate_supported_but_still_conditional(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("constructively forces", candidate["content"])
        self.assertIn("exact multiplicity-ratio certificates", candidate["why_this_is_more_proof_like"])
        self.assertIn("not derived", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
