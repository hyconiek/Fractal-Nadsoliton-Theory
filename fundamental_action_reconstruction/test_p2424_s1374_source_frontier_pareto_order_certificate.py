#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2424_s1374_source_frontier_pareto_order_certificate.py"
OUT = ROOT / "generated" / "p2424_s1374_source_frontier_pareto_order_certificate.json"
MD = ROOT / "generated" / "p2424_s1374_source_frontier_pareto_order_certificate.md"
P2423 = ROOT / "generated" / "p2423_s1373_admissible_repair_order_poset_certificate.json"


class P2424SourceFrontierParetoOrderCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2423.exists():
            subprocess.run([sys.executable, str(ROOT / "p2423_s1373_admissible_repair_order_poset_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["source_frontier_pareto_order_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2424")
        self.assertEqual(self.payload["stage_id"], "S1374")
        self.assertEqual(self.payload["status"], "PASS_SOURCE_FRONTIER_PARETO_ORDER_NO_UNIQUE_SOURCE_GATE_NO_DISCHARGE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_pareto_counts(self) -> None:
        self.assertEqual(self.theorem["source_frontier_gate_count"], 3)
        self.assertEqual(self.theorem["admissible_order_count"], 6)
        self.assertEqual(self.theorem["pareto_order_count"], 4)
        self.assertEqual(self.theorem["dominated_order_count"], 2)
        self.assertEqual(self.theorem["class_counts"], {"bridge_first_pareto": 2, "mixed_split_dominated": 2, "selector_pair_first_pareto": 2})
        self.assertEqual(self.theorem["pareto_class_counts"], {"bridge_first_pareto": 2, "selector_pair_first_pareto": 2})

    def test_pareto_vectors_and_no_unique_first_gate(self) -> None:
        self.assertEqual(self.theorem["objective_vector_counts"], {"[1, 3]": 2, "[2, 3]": 2, "[3, 2]": 2})
        self.assertEqual(self.theorem["pareto_objective_vector_counts"], {"[1, 3]": 2, "[3, 2]": 2})
        self.assertFalse(self.theorem["unique_pareto_first_gate_selected"])

    def test_inheritance_and_hard_limits(self) -> None:
        self.assertTrue(self.theorem["p2423_admissible_order_count_inherited"])
        self.assertTrue(self.theorem["p2423_role_transfer_step_distribution_inherited"])
        self.assertTrue(self.theorem["p2423_ltotal_step_distribution_inherited"])
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["chi11_source_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])

    def test_gatekeepers_and_docs(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("Pareto", MD.read_text(encoding="utf-8"))
        self.assertIn("P2424/S1374", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2424/S1374", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
