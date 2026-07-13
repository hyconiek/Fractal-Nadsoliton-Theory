import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3144_s2094_upsilon_sel_unit_measure_obstruction.py"
OUT = ROOT / "generated" / "p3144_s2094_upsilon_sel_unit_measure_obstruction.json"
MD = ROOT / "generated" / "p3144_s2094_upsilon_sel_unit_measure_obstruction.md"


class P3144UpsilonSelUnitMeasureObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_orbit_counts(self):
        self.assertEqual(self.payload["status"], "P3144_UPSILON_SEL_UNIT_MEASURE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3143"])
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["point_count"], 24)
        self.assertEqual(counts["generator_count"], 3)
        self.assertEqual(counts["orbit_count"], 1)
        self.assertEqual(counts["largest_orbit_size"], 24)
        self.assertEqual(counts["candidate_measures_tested"], 3)
        self.assertEqual(counts["accepted_unit_selector_measures"], 0)

    def test_candidate_measure_obstruction(self):
        rows = {row["candidate_id"]: row for row in self.payload["candidate_measure_rows"]}
        self.assertTrue(rows["U1_uniform_full_orbit"]["unit_normalized"])
        self.assertTrue(rows["U1_uniform_full_orbit"]["full_symmetry_invariant"])
        self.assertFalse(rows["U1_uniform_full_orbit"]["localized_on_selected_branch"])
        self.assertTrue(rows["U2_delta_selected_branch"]["localized_on_selected_branch"])
        self.assertTrue(rows["U2_delta_selected_branch"]["imports_selector_axiom"])
        self.assertFalse(rows["U2_delta_selected_branch"]["full_symmetry_invariant"])
        self.assertTrue(all(not row["accepted_as_upsilon_sel_unit_measure"] for row in rows.values()))

    def test_decision_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertIn("uniform full-orbit measure", decision["bounded_result"])
        self.assertIn("separate conditional", decision["axiomatic_route_recommendation"])
        self.assertIn("Pivot away", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_updated(self):
        self.assertIn("P3144/S2094", MD.read_text(encoding="utf-8"))
        self.assertIn("P3144/S2094", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3144/S2094", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3144", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
