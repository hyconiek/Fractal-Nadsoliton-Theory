#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2420_s1370_bridge_selector_nonclosure_reason_matrix_certificate.py"
OUT = ROOT / "generated" / "p2420_s1370_bridge_selector_nonclosure_reason_matrix_certificate.json"
MD = ROOT / "generated" / "p2420_s1370_bridge_selector_nonclosure_reason_matrix_certificate.md"
P2419 = ROOT / "generated" / "p2419_s1369_chi11_phase_selector_coupling_cut_certificate.json"


class P2420BridgeSelectorNonclosureReasonMatrixCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2419.exists():
            subprocess.run([sys.executable, str(ROOT / "p2419_s1369_chi11_phase_selector_coupling_cut_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["bridge_selector_nonclosure_reason_matrix_certificate"]
        cls.theorem = cls.cert["theorem_export"]
        cls.finite = cls.cert["finite_witness_certificate"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2420")
        self.assertEqual(self.payload["stage_id"], "S1370")
        self.assertEqual(self.payload["status"], "PASS_BRIDGE_SELECTOR_MECHANISM_NONCLOSURE_REASON_MATRIX_NO_TOE_CLOSURE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_closure_truth_table_counts(self) -> None:
        self.assertEqual(self.theorem["closure_gate_count"], 7)
        self.assertEqual(self.theorem["total_closure_assignment_count"], 128)
        self.assertEqual(self.theorem["current_true_gates"], ["apd_value_bridge_witness", "chi11_phase_selector_cut_mechanism"])
        self.assertEqual(self.theorem["current_mask"], 5)
        self.assertEqual(self.theorem["full_mask"], 127)
        self.assertEqual(self.theorem["minimal_toe_mask_count"], 1)
        self.assertEqual(self.theorem["minimal_toe_masks"], [127])
        self.assertEqual(self.theorem["toe_true_mask_count"], 1)

    def test_bridge_selector_subcube_still_nonclosure(self) -> None:
        self.assertEqual(self.theorem["apd_plus_selector_mechanism_subcube_count"], 32)
        self.assertEqual(self.theorem["apd_plus_selector_mechanism_subcube_toe_count"], 1)
        self.assertEqual(self.theorem["apd_plus_selector_mechanism_subcube_failure_count"], 31)
        self.assertEqual(self.theorem["current_to_toe_minimum_repair_distance"], 5)
        self.assertEqual(self.theorem["single_flip_from_current_closes_toe_count"], 0)
        self.assertFalse(self.theorem["current_full_bridge_ready"])
        self.assertFalse(self.theorem["current_selector_closed"])
        self.assertFalse(self.theorem["current_toe_ready"])

    def test_reasons_and_inherited_statuses(self) -> None:
        self.assertEqual(self.theorem["reason_count"], 5)
        self.assertIn("R1_VALUE_BRIDGE_IS_NOT_SOURCE_BRIDGE", self.theorem["why_bridge_plus_selector_mechanism_does_not_close"])
        self.assertIn("R2_SELECTOR_MECHANISM_IS_NOT_SELECTOR_SOURCE", self.theorem["why_bridge_plus_selector_mechanism_does_not_close"])
        self.assertIn("R3_QW2191_IS_STILL_AN_INDEPENDENT_SELECTOR_OBSTRUCTION", self.theorem["why_bridge_plus_selector_mechanism_does_not_close"])
        self.assertTrue(self.theorem["apd_value_bridge_witness_inherited"])
        self.assertTrue(self.theorem["source_discharge_zero_inherited"])
        self.assertTrue(self.theorem["chi11_cut_mechanism_inherited"])

    def test_hard_limits_and_docs(self) -> None:
        self.assertFalse(self.theorem["chi11_source_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("bridge-selector nonclosure", MD.read_text(encoding="utf-8"))
        self.assertIn("P2420/S1370 bridge-selector", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2420/S1370 bridge-selector", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
