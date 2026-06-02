#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2429_s1379_repair_derivative_nearest_miss_witness_table_certificate.py"
OUT = ROOT / "generated" / "p2429_s1379_repair_derivative_nearest_miss_witness_table_certificate.json"
MD = ROOT / "generated" / "p2429_s1379_repair_derivative_nearest_miss_witness_table_certificate.md"
P2428 = ROOT / "generated" / "p2428_s1378_repair_readiness_anf_derivative_certificate.json"


class P2429RepairDerivativeNearestMissWitnessTableCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2428.exists():
            subprocess.run([sys.executable, str(ROOT / "p2428_s1378_repair_readiness_anf_derivative_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["repair_derivative_nearest_miss_witness_table_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2429")
        self.assertEqual(self.payload["stage_id"], "S1379")
        self.assertEqual(self.payload["status"], "PASS_REPAIR_DERIVATIVE_WITNESS_TABLE_NO_GATE_DISCHARGE_NO_CLOSURE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_witness_counts(self) -> None:
        self.assertEqual(self.theorem["truth_table_size"], 32)
        self.assertEqual(self.theorem["witness_row_count"], 69)
        self.assertEqual(
            self.theorem["witness_count_by_target"],
            {
                "bridge_source_ready": 16,
                "role_bearing_ltotal_ready": 16,
                "role_transfer_ready": 16,
                "selector_source_ready": 16,
                "toe_ready": 5,
            },
        )
        self.assertEqual(self.theorem["target_gate_pair_with_witness_count"], 10)

    def test_nearest_miss_structure(self) -> None:
        self.assertEqual(self.theorem["toe_nearest_miss_row_count"], 5)
        self.assertTrue(self.theorem["toe_nearest_miss_one_per_gate"])
        self.assertEqual(self.theorem["toe_nearest_miss_missing_before_distribution"], {"1": 5})
        self.assertTrue(self.theorem["selector_derivative_witnesses_only_chi11_qw2191"])
        self.assertTrue(self.theorem["role_bearing_ltotal_witnesses_only_ltotal_export"])

    def test_inherited_hard_limits_and_docs(self) -> None:
        self.assertTrue(self.theorem["p2428_derivative_edge_counts_inherited"])
        self.assertTrue(self.theorem["p2428_toe_derivative_edges_inherited"])
        self.assertTrue(self.theorem["p2428_selector_pair_inherited"])
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["chi11_source_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("nearest-miss witness", MD.read_text(encoding="utf-8"))
        self.assertIn("P2429/S1379", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2429/S1379", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
