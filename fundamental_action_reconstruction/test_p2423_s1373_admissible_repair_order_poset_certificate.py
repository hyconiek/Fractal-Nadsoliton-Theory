#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2423_s1373_admissible_repair_order_poset_certificate.py"
OUT = ROOT / "generated" / "p2423_s1373_admissible_repair_order_poset_certificate.json"
MD = ROOT / "generated" / "p2423_s1373_admissible_repair_order_poset_certificate.md"
P2422 = ROOT / "generated" / "p2422_s1372_current_missing_gate_repair_subcube_certificate.json"


class P2423AdmissibleRepairOrderPosetCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2422.exists():
            subprocess.run([sys.executable, str(ROOT / "p2422_s1372_current_missing_gate_repair_subcube_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["admissible_repair_order_poset_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2423")
        self.assertEqual(self.payload["stage_id"], "S1373")
        self.assertEqual(self.payload["status"], "PASS_ADMISSIBLE_REPAIR_ORDER_POSET_NO_GATE_DISCHARGE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_order_counts(self) -> None:
        self.assertEqual(self.theorem["missing_gate_count"], 5)
        self.assertEqual(self.theorem["precedence_edge_count"], 4)
        self.assertEqual(self.theorem["raw_order_count"], 120)
        self.assertEqual(self.theorem["admissible_order_count"], 6)
        self.assertEqual(self.theorem["rejected_order_count"], 114)

    def test_step_distributions(self) -> None:
        self.assertEqual(self.theorem["role_transfer_ready_step_distribution"], {"4": 6})
        self.assertEqual(self.theorem["role_bearing_ltotal_ready_step_distribution"], {"5": 6})
        self.assertEqual(self.theorem["toe_ready_step_distribution"], {"5": 6})
        self.assertEqual(self.theorem["selector_source_ready_step_distribution"], {"2": 2, "3": 4})
        self.assertTrue(self.theorem["all_edges_necessary"])

    def test_inheritance_and_hard_limits(self) -> None:
        self.assertTrue(self.theorem["p2422_repair_subcube_count_inherited"])
        self.assertTrue(self.theorem["p2422_toe_ready_repair_count_inherited"])
        self.assertTrue(self.theorem["p2422_selector_singleton_unlock_count_inherited"])
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["chi11_source_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("repair-order poset", MD.read_text(encoding="utf-8"))
        self.assertIn("P2423/S1373", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2423/S1373", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
