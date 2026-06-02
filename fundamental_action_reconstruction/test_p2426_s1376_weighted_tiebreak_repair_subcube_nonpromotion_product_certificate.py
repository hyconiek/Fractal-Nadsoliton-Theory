#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2426_s1376_weighted_tiebreak_repair_subcube_nonpromotion_product_certificate.py"
OUT = ROOT / "generated" / "p2426_s1376_weighted_tiebreak_repair_subcube_nonpromotion_product_certificate.json"
MD = ROOT / "generated" / "p2426_s1376_weighted_tiebreak_repair_subcube_nonpromotion_product_certificate.md"
P2425 = ROOT / "generated" / "p2425_s1375_source_frontier_weighted_tiebreak_premise_certificate.json"


class P2426WeightedTiebreakRepairSubcubeNonpromotionProductCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2425.exists():
            subprocess.run([sys.executable, str(ROOT / "p2425_s1375_source_frontier_weighted_tiebreak_premise_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["weighted_tiebreak_repair_subcube_nonpromotion_product_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2426")
        self.assertEqual(self.payload["stage_id"], "S1376")
        self.assertEqual(self.payload["status"], "PASS_WEIGHTED_REPAIR_PRODUCT_NO_GATE_DISCHARGE_NO_CLOSURE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_product_counts(self) -> None:
        self.assertEqual(self.theorem["weight_assignment_count"], 144)
        self.assertEqual(self.theorem["repair_subcube_assignment_count"], 32)
        self.assertEqual(self.theorem["product_assignment_count"], 4608)
        self.assertEqual(self.theorem["toe_ready_product_count"], 144)
        self.assertEqual(self.theorem["proper_repair_failure_product_count"], 4464)
        self.assertEqual(self.theorem["empty_repair_row_count"], 144)
        self.assertEqual(self.theorem["empty_repair_remaining_missing_distribution"], {"5": 144})

    def test_weight_counts_and_inheritance(self) -> None:
        self.assertEqual(
            self.theorem["weight_winner_product_counts"],
            {"bridge_first_pareto": 3456, "bridge_selector_tie": 192, "selector_pair_first_pareto": 960},
        )
        self.assertEqual(
            self.theorem["toe_ready_count_by_weight_winner"],
            {"bridge_first_pareto": 108, "bridge_selector_tie": 6, "selector_pair_first_pareto": 30},
        )
        self.assertTrue(self.theorem["p2425_winner_counts_inherited"])
        self.assertTrue(self.theorem["p2425_no_internal_weight_premise_inherited"])
        self.assertTrue(self.theorem["p2422_repair_subcube_count_inherited"])
        self.assertTrue(self.theorem["p2422_toe_ready_repair_count_inherited"])
        self.assertFalse(self.theorem["weighted_choice_reduces_missing_gate_count"])

    def test_hard_limits_and_docs(self) -> None:
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["chi11_source_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("nonpromotion product", MD.read_text(encoding="utf-8"))
        self.assertIn("P2426/S1376", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2426/S1376", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
