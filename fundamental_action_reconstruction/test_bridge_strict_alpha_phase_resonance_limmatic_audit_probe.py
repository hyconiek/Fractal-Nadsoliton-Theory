#!/usr/bin/env python3
"""Tests for the strict-alpha phase resonance limmatic audit probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_phase_resonance_limmatic_audit_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_phase_resonance_limmatic_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_phase_resonance_limmatic_audit_report.md"


class BridgeStrictAlphaPhaseResonanceLimmaticAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(self.payload["result_kind"], "SCRATCH_STRICT_ALPHA_PHASE_RESONANCE_LIMMATIC_AUDIT_PROBE__NOT_A_THEOREM")
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("The 12-fifth closure residual is 3^12/2^19, not 256/243.", self.payload["hard_limits"])
        self.assertIn("No Aut(Z_12)-invariant generator/orientation fixing is claimed; phase labelling remains premise-based under the N462 boundary.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_limma_identities_are_true(self) -> None:
        limma = self.payload["limma_identity_check"]
        self.assertEqual(limma["q_radical_power_5"], "256/243")
        self.assertTrue(limma["limma_is_2_8_over_3_5"])
        self.assertEqual(limma["limma_from_fourth_div_two_whole_tones"], "256/243")
        self.assertTrue(limma["limma_from_fourth_identity_ok"])
        self.assertEqual(limma["limma_times_apotome_equals_whole_tone"], "9/8")
        self.assertTrue(limma["limma_apotome_identity_ok"])
        self.assertLess(abs(limma["eta_residual_vs_9_5"]), 1e-12)

    def test_twelve_fifth_closure_claim_is_rejected(self) -> None:
        closure = self.payload["closure_residual_discriminator"]
        self.assertEqual(closure["actual_12_fifth_closure_residual"], "531441/524288")
        self.assertEqual(closure["actual_12_fifth_closure_name"], "Pythagorean comma")
        self.assertEqual(closure["limma_fraction"], "256/243")
        self.assertFalse(closure["limma_equals_12_fifth_closure_residual"])
        self.assertIn("REJECT_12_FIFTH_CLOSURE_EQUALS_LIMMA", closure["verdict"])

    def test_z12_phase_labelling_is_conditional(self) -> None:
        z12 = self.payload["twelve_node_resonance_audit"]
        self.assertEqual(z12["z12_fifth_step"], 7)
        self.assertTrue(z12["fifth_step_is_generator_mod_12"])
        self.assertTrue(z12["fifth_cycle_visits_all_nodes"])
        self.assertEqual(sorted(z12["fifth_cycle_mod_12"]), list(range(12)))
        self.assertTrue(z12["fifth_step_equivalent_to_generator_only_after_generator_gauge_choice"])
        phase = self.payload["phase_labelled_selector_interpretation"]
        self.assertTrue(phase["supported_conditionally"])
        self.assertIn("N462 boundary", phase["strict_core_blocker"])

    def test_candidate_supported_but_not_closure(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("arithmetically real", candidate["content"])
        self.assertIn("false 12-fifth closure", candidate["why_this_is_more_proof_like"])
        self.assertIn("does not derive", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
