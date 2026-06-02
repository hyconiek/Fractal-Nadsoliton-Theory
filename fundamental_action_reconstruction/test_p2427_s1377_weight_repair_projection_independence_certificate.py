#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2427_s1377_weight_repair_projection_independence_certificate.py"
OUT = ROOT / "generated" / "p2427_s1377_weight_repair_projection_independence_certificate.json"
MD = ROOT / "generated" / "p2427_s1377_weight_repair_projection_independence_certificate.md"
P2426 = ROOT / "generated" / "p2426_s1376_weighted_tiebreak_repair_subcube_nonpromotion_product_certificate.json"


class P2427WeightRepairProjectionIndependenceCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2426.exists():
            subprocess.run(
                [sys.executable, str(ROOT / "p2426_s1376_weighted_tiebreak_repair_subcube_nonpromotion_product_certificate.py")],
                check=True,
            )
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["weight_repair_projection_independence_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2427")
        self.assertEqual(self.payload["stage_id"], "S1377")
        self.assertEqual(
            self.payload["status"],
            "PASS_WEIGHT_REPAIR_PROJECTION_INDEPENDENCE_NO_GATE_DISCHARGE_NO_CLOSURE",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_projection_product_counts(self) -> None:
        self.assertEqual(self.theorem["weight_assignment_count"], 144)
        self.assertEqual(self.theorem["repair_assignment_count"], 32)
        self.assertEqual(self.theorem["product_assignment_count"], 4608)
        self.assertEqual(
            self.theorem["weight_counts"],
            {"bridge_first_pareto": 108, "bridge_selector_tie": 6, "selector_pair_first_pareto": 30},
        )
        self.assertEqual(
            self.theorem["repair_missing_count_distribution"],
            {"0": 1, "1": 5, "2": 10, "3": 10, "4": 5, "5": 1},
        )

    def test_readiness_independence(self) -> None:
        self.assertEqual(
            self.theorem["repair_readiness_true_counts"],
            {
                "bridge_source_ready": 16,
                "role_bearing_ltotal_ready": 16,
                "role_transfer_ready": 16,
                "selector_source_ready": 8,
                "toe_ready": 1,
            },
        )
        self.assertTrue(self.theorem["all_readiness_tables_factor_by_weight_x_repair"])
        self.assertTrue(self.theorem["conditional_readiness_distributions_identical_across_weight_winners"])
        self.assertTrue(self.theorem["conditional_missing_count_distribution_identical_across_weight_winners"])
        self.assertEqual(self.theorem["weight_side_support_for_repair_readiness_predicates"], [])
        self.assertEqual(self.theorem["repair_side_support_for_weight_winner_predicate"], [])

    def test_inherited_hard_limits_and_docs(self) -> None:
        self.assertTrue(self.theorem["p2426_product_assignment_count_inherited"])
        self.assertTrue(self.theorem["p2426_weight_winner_counts_inherited"])
        self.assertTrue(self.theorem["p2426_toe_ready_product_count_inherited"])
        self.assertTrue(self.theorem["p2426_weighted_choice_reduces_missing_gate_count_inherited"])
        self.assertFalse(self.theorem["weighted_choice_discharges_repair_gate"])
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["chi11_source_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("projection independence", MD.read_text(encoding="utf-8"))
        self.assertIn("P2427/S1377", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2427/S1377", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
