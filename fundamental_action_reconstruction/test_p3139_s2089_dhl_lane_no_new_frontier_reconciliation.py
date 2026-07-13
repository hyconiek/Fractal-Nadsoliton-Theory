import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3139_s2089_dhl_lane_no_new_frontier_reconciliation.py"
OUT = ROOT / "generated" / "p3139_s2089_dhl_lane_no_new_frontier_reconciliation.json"
MD = ROOT / "generated" / "p3139_s2089_dhl_lane_no_new_frontier_reconciliation.md"


class P3139DHLLaneNoNewFrontierReconciliationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P3139_DHL_LANE_NO_NEW_LIVE_FRONTIER_RECONCILIATION")
        self.assertTrue(all(self.payload["input_hashes"].values()))
        self.assertGreaterEqual(self.payload["repo_grep_summary"]["D_HL"], 1)
        self.assertGreaterEqual(self.payload["repo_grep_summary"]["QW-2191"], 1)

    def test_obligation_matrix(self):
        obligations = {row["obligation"]: row["satisfied_by_current_lane"] for row in self.payload["source_obligations"]}
        self.assertTrue(obligations["O1_explicit_object"])
        self.assertTrue(obligations["O2_nonzero_inversion_odd_value"])
        self.assertFalse(obligations["O3_absolute_support_origin_after_Z12_quotient"])
        self.assertFalse(obligations["O4_unpaired_polarity_lambda"])
        self.assertFalse(obligations["O5_import_free_strict_source_law"])
        self.assertFalse(obligations["O6_variational_unit_coupling"])

    def test_theorem_and_decision(self):
        theorem = self.payload["theorem"]
        self.assertEqual(theorem["finite_support"]["P3138_profiles_tested"], 120)
        self.assertEqual(theorem["finite_support"]["P3138_translation_covariant_receiver_rows"], 120)
        self.assertEqual(theorem["finite_support"]["P3138_accepted_import_free_joint_sources"], 0)
        self.assertIn("receiver-to-source", self.payload["constructed_object"]["name"])
        self.assertEqual(self.payload["decision"]["accepted_import_free_source_rows"], 0)
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_docs_updated(self):
        self.assertIn("P3139/S2089", MD.read_text(encoding="utf-8"))
        self.assertIn("P3139/S2089", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3139/S2089", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3139", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
