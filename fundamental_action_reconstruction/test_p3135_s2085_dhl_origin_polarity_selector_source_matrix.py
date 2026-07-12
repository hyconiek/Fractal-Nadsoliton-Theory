import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3135_s2085_dhl_origin_polarity_selector_source_matrix.py"
OUT = ROOT / "generated" / "p3135_s2085_dhl_origin_polarity_selector_source_matrix.json"
MD = ROOT / "generated" / "p3135_s2085_dhl_origin_polarity_selector_source_matrix.md"


class P3135DHLOriginPolaritySelectorSourceMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_backscan(self):
        self.assertEqual(self.payload["status"], "P3135_DHL_JOINT_ORIGIN_POLARITY_SELECTOR_SOURCE_MATRIX_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3134"])
        self.assertGreaterEqual(len(self.payload["repo_grep_backscan_summary"]), 3)

    def test_pair_orbit_obstruction(self):
        action = self.payload["finite_group_action"]
        self.assertEqual(action["pair_count"], 24)
        self.assertEqual(action["orbit_count"], 1)
        self.assertEqual(len(action["orbits"][0]["members"]), 24)
        self.assertEqual(action["invariant_unique_pair_selectors"], 0)

    def test_source_matrix_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["candidate_sources_tested"], 10)
        self.assertEqual(cert["joint_pair_orbits"], 1)
        self.assertEqual(cert["largest_pair_orbit"], 24)
        self.assertEqual(cert["accepted_joint_r_lambda_sources"], 0)
        self.assertEqual(cert["sources_passing_all_three_gates"], 0)
        rows = self.payload["source_acceptance_rows"]
        self.assertTrue(any(row["candidate_source"] == "P2718_chiral_bispectrum_sign" and row["selects_lambda"] for row in rows))
        self.assertTrue(any(row["candidate_source"] == "P3130_Theta_TO_translation_quotient" and row["import_free_current_artifacts"] for row in rows))
        self.assertTrue(all(not row["accepted_joint_r_lambda_source"] for row in rows))

    def test_decision_and_docs(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["repo_backscan_used"])
        self.assertTrue(decision["positive_scoped_flags"]["finite_orbit_obstruction_proved"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("J_DHL", decision["next_honest_step"])
        self.assertIn("P3135/S2085", MD.read_text(encoding="utf-8"))
        self.assertIn("P3135/S2085", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3135/S2085", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3135", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
