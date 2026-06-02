#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2422_s1372_current_missing_gate_repair_subcube_certificate.py"
OUT = ROOT / "generated" / "p2422_s1372_current_missing_gate_repair_subcube_certificate.json"
MD = ROOT / "generated" / "p2422_s1372_current_missing_gate_repair_subcube_certificate.md"
P2421 = ROOT / "generated" / "p2421_s1371_bridge_selector_closure_prime_implicant_failure_cover_certificate.json"


class P2422CurrentMissingGateRepairSubcubeCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2421.exists():
            subprocess.run([sys.executable, str(ROOT / "p2421_s1371_bridge_selector_closure_prime_implicant_failure_cover_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["current_missing_gate_repair_subcube_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2422")
        self.assertEqual(self.payload["stage_id"], "S1372")
        self.assertEqual(self.payload["status"], "PASS_CURRENT_REPAIR_SUBCUBE_PARTIAL_UNLOCKS_NO_GATE_DISCHARGE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_repair_subcube_counts(self) -> None:
        self.assertEqual(self.theorem["missing_gate_count"], 5)
        self.assertEqual(self.theorem["repair_subcube_assignment_count"], 32)
        self.assertEqual(self.theorem["proper_repair_failure_count"], 31)
        self.assertEqual(self.theorem["toe_ready_repair_count"], 1)
        self.assertEqual(
            self.theorem["toe_ready_added_gates"],
            [
                "source_obligation_discharge",
                "chi11_source_export",
                "qw2191_selector_discharge",
                "role_transfer_audit_license",
                "role_bearing_ltotal_export",
            ],
        )

    def test_minimal_partial_unlocks(self) -> None:
        self.assertEqual(self.theorem["bridge_source_minimal_unlock_sets"][0]["added_gates"], ["source_obligation_discharge"])
        self.assertEqual(
            self.theorem["selector_source_minimal_unlock_sets"][0]["added_gates"],
            ["chi11_source_export", "qw2191_selector_discharge"],
        )
        self.assertEqual(self.theorem["role_transfer_minimal_unlock_sets"][0]["added_gates"], ["role_transfer_audit_license"])
        self.assertEqual(self.theorem["role_bearing_ltotal_minimal_unlock_sets"][0]["added_gates"], ["role_bearing_ltotal_export"])
        self.assertEqual(self.theorem["selector_singleton_unlock_count"], 0)
        self.assertEqual(self.theorem["singleton_non_toe_unlock_count"], 3)

    def test_inheritance_and_hard_limits(self) -> None:
        self.assertTrue(self.theorem["p2421_current_missing_gates_inherited"])
        self.assertTrue(self.theorem["p2421_repair_distance_inherited"])
        self.assertTrue(self.theorem["p2421_unique_prime_implicant_inherited"])
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["chi11_source_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("repair subcube", MD.read_text(encoding="utf-8"))
        self.assertIn("P2422/S1372", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2422/S1372", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
