import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2773_s1723_shannon_entropy_geometry_forcing_obstruction.py"
OUT = ROOT / "generated" / "p2773_s1723_shannon_entropy_geometry_forcing_obstruction.json"
MD = ROOT / "generated" / "p2773_s1723_shannon_entropy_geometry_forcing_obstruction.md"


class P2773ShannonEntropyGeometryForcingObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2773_SHANNON_ENTROPY_GEOMETRY_FORCING_OBSTRUCTION_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2772"], "P2772_SELF_LEARNING_KERNEL_UPDATE_LAW_STATIONARITY_WITNESS_BOUNDED_NO_GO_NO_CLOSURE")
        self.assertEqual(self.payload["audited_question"], "Does Shannon entropy alone force the nadsoliton geometry?")

    def test_entropy_matches_alpha_geo_but_geometries_differ(self):
        witness = self.payload["entropy_geometry_witness"]
        self.assertEqual(witness["support_size"], 16)
        self.assertTrue(witness["entropy_matches_alpha_geo"])
        self.assertEqual(witness["geometry_row_count"], 4)
        self.assertGreater(witness["distinct_distance_histogram_count"], 1)
        self.assertFalse(witness["entropy_forces_unique_geometry_on_this_class"])
        geometries = {row["geometry"] for row in witness["geometry_rows"]}
        self.assertEqual(geometries, {"complete", "cycle", "path", "star"})
        diameters = {row["diameter"] for row in witness["geometry_rows"]}
        self.assertGreater(len(diameters), 1)

    def test_acceptance_blocks_entropy_geometry_forcing(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["facts"]["entropy_matches_alpha_geo"])
        self.assertTrue(acceptance["facts"]["multiple_inequivalent_geometries_at_same_entropy"])
        self.assertFalse(acceptance["facts"]["entropy_forces_unique_geometry_on_test_class"])
        self.assertFalse(acceptance["accepted_as_shannon_entropy_geometry_forcing_theorem"])
        self.assertFalse(acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(acceptance["accepted_as_kernel_full_expression_theorem"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("entropy_forces_unique_geometry_on_test_class", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("entropy-plus", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2773", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2773/S1723", MD.read_text(encoding="utf-8"))
        self.assertIn("P2773/S1723", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2773/S1723", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2773", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
