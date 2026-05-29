#!/usr/bin/env python3
"""Tests for the strict-alpha distance-5 band-pass support selector certificate probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_bandpass_d5_support_selector_certificate_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_bandpass_d5_support_selector_certificate_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_bandpass_d5_support_selector_certificate_report.md"


class BridgeStrictAlphaBandpassD5SupportSelectorCertificateProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(
            self.payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_BANDPASS_D5_SUPPORT_SELECTOR_CERTIFICATE_PROBE__NOT_A_THEOREM",
        )
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("The distance-5 band-pass support selector is a conditional model-internal support premise, not a derived strict nadsoliton action theorem.", self.payload["hard_limits"])
        self.assertIn("The fifth-step orbit is selected only up to translate; no unique source phase, orientation, or assignment theorem is exported.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_previous_polarity_and_target_replay(self) -> None:
        replay = self.payload["previous_polarity_replay"]
        self.assertIn("SUPPORT_DISTANCE_HISTOGRAM_POLARITY", replay["result_kind"])
        self.assertEqual(replay["plain_dirichlet_gap_formula"], "E_ind(fifth)-E_ind(contiguous)=8*(w_1-w_5)")
        self.assertTrue(replay["plain_dirichlet_prefers_contiguous"])
        self.assertTrue(replay["plain_dirichlet_makes_fifth_maximum"])

        target = self.payload["target_identity_replay"]
        self.assertEqual(target["q_power_product"], "256/243")
        self.assertLess(abs(target["eta_residual_vs_9_5"]), 1e-12)
        self.assertEqual(target["balanced_ledger"], [2, 2, 2, 1, 1])

    def test_distance_5_cycle_theorem(self) -> None:
        theorem = self.payload["distance_5_cycle_theorem"]
        self.assertEqual(theorem["gcd_5_12"], 1)
        self.assertTrue(theorem["is_single_12_cycle"])
        self.assertEqual(theorem["support_size"], 5)
        self.assertEqual(theorem["max_induced_edges_bound"], 4)
        self.assertIn("proper k-vertex subset of a cycle", theorem["reason"])

    def test_bandpass_scan_selects_fifth_step_orbit(self) -> None:
        scan = self.payload["bandpass_selector_scan"]
        self.assertEqual(scan["support_count"], 792)
        self.assertEqual(scan["max_h5"], 4)
        self.assertEqual(scan["next_best_h5"], 3)
        self.assertEqual(scan["h5_gap_to_next_best"], 1)
        self.assertEqual(scan["min_bandpass_energy_minus_h5"], -4)
        self.assertEqual(scan["maximizer_count"], 12)
        self.assertTrue(scan["maximizers_equal_fifth_step_orbit"])
        self.assertEqual(scan["fifth_step_support_ordered"], [0, 7, 2, 9, 4])
        self.assertEqual(scan["fifth_step_support_canonical"], [0, 2, 4, 7, 9])
        self.assertEqual(scan["fifth_step_histogram"], [0, 3, 2, 1, 4, 0])
        self.assertEqual(scan["contiguous_histogram"], [4, 3, 2, 1, 0, 0])
        self.assertEqual(scan["fifth_step_h5"], 4)
        self.assertEqual(scan["contiguous_h5"], 0)
        self.assertIn([0, 2, 4, 7, 9], [row["support"] for row in scan["maximizers"]])

    def test_conditional_schema_and_candidate_status(self) -> None:
        schema = self.payload["conditional_selector_schema"]
        self.assertIn("maximizes h_5", schema["support_premise"])
        self.assertIn("balanced ledger", schema["ledger_premise"])
        self.assertIn("fifth-step orbit", schema["conditional_output"])
        self.assertTrue(any("choose a unique translate/source phase" in item for item in schema["remaining_selector_obligations"]))
        self.assertTrue(any("discharge QW-2191" in item for item in schema["remaining_selector_obligations"]))

        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("maximizing distance-5 internal pairs", candidate["content"])
        self.assertIn("distance-5 graph being a 12-cycle", candidate["why_this_is_more_proof_like"])
        self.assertIn("does not derive the required distance-5 band-pass action", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
