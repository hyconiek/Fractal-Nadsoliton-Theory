#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2428_s1378_repair_readiness_anf_derivative_certificate.py"
OUT = ROOT / "generated" / "p2428_s1378_repair_readiness_anf_derivative_certificate.json"
MD = ROOT / "generated" / "p2428_s1378_repair_readiness_anf_derivative_certificate.md"
P2427 = ROOT / "generated" / "p2427_s1377_weight_repair_projection_independence_certificate.json"


class P2428RepairReadinessANFDerivativeCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2427.exists():
            subprocess.run([sys.executable, str(ROOT / "p2427_s1377_weight_repair_projection_independence_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["repair_readiness_anf_derivative_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2428")
        self.assertEqual(self.payload["stage_id"], "S1378")
        self.assertEqual(self.payload["status"], "PASS_REPAIR_READINESS_ANF_DERIVATIVE_NO_GATE_DISCHARGE_NO_CLOSURE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_anf_prime_implicants(self) -> None:
        self.assertEqual(self.theorem["gate_count"], 5)
        self.assertEqual(self.theorem["truth_table_size"], 32)
        self.assertTrue(self.theorem["single_prime_implicant_per_readiness_predicate"])
        self.assertTrue(self.theorem["selector_requires_chi11_and_qw2191_pair"])
        self.assertTrue(self.theorem["toe_prime_implicant_is_all_five_gates"])
        self.assertEqual(self.theorem["anf_degrees"]["toe_ready"], 5)
        self.assertEqual(self.theorem["anf_degrees"]["selector_source_ready"], 2)

    def test_derivative_edge_counts(self) -> None:
        self.assertEqual(
            self.theorem["derivative_edge_counts"]["bridge_source_ready"],
            {
                "source_obligation_discharge": 16,
                "chi11_source_export": 0,
                "qw2191_selector_discharge": 0,
                "role_transfer_audit_license": 0,
                "role_bearing_ltotal_export": 0,
            },
        )
        self.assertEqual(
            self.theorem["derivative_edge_counts"]["selector_source_ready"],
            {
                "source_obligation_discharge": 0,
                "chi11_source_export": 8,
                "qw2191_selector_discharge": 8,
                "role_transfer_audit_license": 0,
                "role_bearing_ltotal_export": 0,
            },
        )
        self.assertTrue(self.theorem["toe_derivative_edges_are_nearest_misses"])

    def test_inherited_hard_limits_and_docs(self) -> None:
        self.assertTrue(self.theorem["p2427_product_assignment_count_inherited"])
        self.assertTrue(self.theorem["p2427_all_tables_factor_inherited"])
        self.assertTrue(self.theorem["p2427_weight_side_support_empty_inherited"])
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["chi11_source_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("ANF derivative", MD.read_text(encoding="utf-8"))
        self.assertIn("P2428/S1378", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2428/S1378", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
