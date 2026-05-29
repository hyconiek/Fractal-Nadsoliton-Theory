#!/usr/bin/env python3
"""Tests for the strict-alpha strict-gate Dirichlet action audit probe."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_strict_gate_dirichlet_action_audit_probe.py"
REPORT = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_strict_gate_dirichlet_action_audit_report.json"
REPORT_MD = ROOT / "fundamental_action_reconstruction" / "scratch" / "bridge_strict_alpha_strict_gate_dirichlet_action_audit_report.md"


class BridgeStrictAlphaStrictGateDirichletActionAuditProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(REPORT.read_text(encoding="utf-8"))

    def test_result_kind_and_guardrails(self) -> None:
        self.assertEqual(
            self.payload["result_kind"],
            "SCRATCH_STRICT_ALPHA_STRICT_GATE_DIRICHLET_ACTION_AUDIT_PROBE__NOT_A_THEOREM",
        )
        self.assertIn("No identity K_legacy_ont == K_strict_gate is exported.", self.payload["hard_limits"])
        self.assertIn("The strict-gate Dirichlet quadratic action is audited as a candidate construction, not derived as a strict nadsoliton action theorem.", self.payload["hard_limits"])
        self.assertIn("The support scan prefers clustered supports over the fifth-step resonance window; no support/phase-placement theorem is exported.", self.payload["hard_limits"])
        self.assertIn("No QW-2191 selector discharge and no ToE closure are claimed.", self.payload["hard_limits"])

    def test_strict_gate_weight_and_laplacian_positivity(self) -> None:
        params = self.payload["strict_gate_parameters"]
        self.assertEqual(params["kernel"], "K_strict_gate(d)=cos(omega*d+phi)/(1+beta*d^eta)")
        self.assertAlmostEqual(params["omega"], 0.18575)
        self.assertAlmostEqual(params["phi"], 0.16250)
        self.assertAlmostEqual(params["eta"], 9.0 / 5.0)
        weights = self.payload["positive_weight_audit"]
        self.assertTrue(weights["all_cyclic_distance_weights_positive"])
        self.assertEqual([row["distance"] for row in weights["distance_weight_rows"]], [1, 2, 3, 4, 5, 6])
        self.assertTrue(all(row["weight"] > 0.0 for row in weights["distance_weight_rows"]))

        spectrum = self.payload["laplacian_psd_audit"]
        self.assertAlmostEqual(spectrum["zero_mode_eigenvalue"], 0.0)
        self.assertTrue(spectrum["all_nonconstant_modes_positive"])
        self.assertAlmostEqual(spectrum["minimum_nonconstant_eigenvalue"], 0.7541211542070796)
        self.assertEqual([row["mode"] for row in spectrum["spectrum_rows"]], list(range(12)))

    def test_target_identity_and_global_scan(self) -> None:
        target = self.payload["target_identity_replay"]
        self.assertEqual(target["q_power_product"], "256/243")
        self.assertLess(abs(target["eta_residual_vs_9_5"]), 1e-12)
        self.assertEqual(target["canonical_ledgers"], [[4, 1, 1, 1, 1], [3, 2, 1, 1, 1], [2, 2, 2, 1, 1]])

        scan = self.payload["canonical_ledger_assignment_scan"]
        self.assertEqual(scan["total_rows_scanned"], 27720)
        self.assertTrue(scan["global_best_is_balanced_ledger"])
        self.assertTrue(scan["global_best_support_is_cyclic_contiguous"])
        self.assertFalse(scan["global_best_support_is_fifth_window"])
        self.assertEqual(scan["global_best_among_canonical_ledgers"]["ledger"], [2, 2, 2, 1, 1])
        self.assertEqual(scan["global_best_among_canonical_ledgers"]["support"], [0, 1, 2, 10, 11])
        self.assertAlmostEqual(scan["global_best_among_canonical_ledgers"]["energy"], 8.066461680488555)

    def test_ledger_summaries_and_support_obstruction(self) -> None:
        scan = self.payload["canonical_ledger_assignment_scan"]
        summaries = {tuple(row["ledger"]): row for row in scan["ledger_summaries"]}
        self.assertEqual(summaries[(4, 1, 1, 1, 1)]["rows_scanned"], 3960)
        self.assertEqual(summaries[(3, 2, 1, 1, 1)]["rows_scanned"], 15840)
        self.assertEqual(summaries[(2, 2, 2, 1, 1)]["rows_scanned"], 7920)
        self.assertGreater(summaries[(4, 1, 1, 1, 1)]["min_energy"], summaries[(3, 2, 1, 1, 1)]["min_energy"])
        self.assertGreater(summaries[(3, 2, 1, 1, 1)]["min_energy"], summaries[(2, 2, 2, 1, 1)]["min_energy"])

        support = self.payload["support_comparison"]
        self.assertEqual(support["fifth_window_support"], [0, 7, 2, 9, 4])
        self.assertEqual(support["contiguous_support"], [0, 1, 2, 3, 4])
        self.assertIn("separate support/phase-placement theorem", support["verdict"])
        fifth = {tuple(row["ledger"]): row for row in support["fifth_window_rows"]}
        contiguous = {tuple(row["ledger"]): row for row in support["contiguous_rows"]}
        self.assertGreater(fifth[(2, 2, 2, 1, 1)]["min_energy"], contiguous[(2, 2, 2, 1, 1)]["min_energy"])
        self.assertAlmostEqual(fifth[(2, 2, 2, 1, 1)]["min_energy"], 18.29375613429685)
        self.assertAlmostEqual(contiguous[(2, 2, 2, 1, 1)]["min_energy"], 8.066461680488556)

    def test_candidate_supported_but_not_closure(self) -> None:
        candidate = self.payload["candidate_interpretation"]
        self.assertTrue(candidate["supported_by_this_probe"])
        self.assertIn("candidate positive Dirichlet quadratic source", candidate["content"])
        self.assertIn("actual strict-gate parameters", candidate["why_this_is_more_proof_like"])
        self.assertIn("support-selection problem", candidate["why_this_is_not_enough"])
        self.assertTrue(REPORT_MD.exists())


if __name__ == "__main__":
    unittest.main()
