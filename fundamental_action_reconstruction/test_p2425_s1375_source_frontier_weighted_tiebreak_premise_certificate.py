#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2425_s1375_source_frontier_weighted_tiebreak_premise_certificate.py"
OUT = ROOT / "generated" / "p2425_s1375_source_frontier_weighted_tiebreak_premise_certificate.json"
MD = ROOT / "generated" / "p2425_s1375_source_frontier_weighted_tiebreak_premise_certificate.md"
P2424 = ROOT / "generated" / "p2424_s1374_source_frontier_pareto_order_certificate.json"


class P2425SourceFrontierWeightedTiebreakPremiseCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2424.exists():
            subprocess.run([sys.executable, str(ROOT / "p2424_s1374_source_frontier_pareto_order_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["source_frontier_weighted_tiebreak_premise_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2425")
        self.assertEqual(self.payload["stage_id"], "S1375")
        self.assertEqual(self.payload["status"], "PASS_WEIGHTED_TIEBREAK_CONDITIONAL_NO_SOURCE_PREMISE_NO_DISCHARGE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_weight_grid_counts_and_conditions(self) -> None:
        self.assertEqual(self.theorem["objective_class_count"], 3)
        self.assertEqual(self.theorem["pareto_objective_class_count"], 2)
        self.assertEqual(self.theorem["weight_grid_max"], 12)
        self.assertEqual(self.theorem["weight_grid_assignment_count"], 144)
        self.assertEqual(
            self.theorem["weight_grid_winner_counts"],
            {"bridge_first_pareto": 108, "bridge_selector_tie": 6, "selector_pair_first_pareto": 30},
        )
        self.assertEqual(self.theorem["bridge_first_condition"], "w_selector < 2*w_bridge")
        self.assertEqual(self.theorem["selector_pair_first_condition"], "w_selector > 2*w_bridge")
        self.assertEqual(self.theorem["tie_boundary"], "w_selector = 2*w_bridge")
        self.assertEqual(self.theorem["dominated_win_count"], 0)

    def test_inheritance_and_no_unique_premise(self) -> None:
        self.assertTrue(self.theorem["p2424_pareto_order_count_inherited"])
        self.assertTrue(self.theorem["p2424_dominated_order_count_inherited"])
        self.assertTrue(self.theorem["p2424_no_unique_first_gate_inherited"])
        self.assertFalse(self.theorem["internal_weight_source_premise_exported"])
        self.assertFalse(self.theorem["unique_first_gate_selected"])

    def test_hard_limits_and_docs(self) -> None:
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["chi11_source_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("weighted tie-break", MD.read_text(encoding="utf-8"))
        self.assertIn("P2425/S1375", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2425/S1375", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
