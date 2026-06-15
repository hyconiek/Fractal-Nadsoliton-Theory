import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2774_s1724_entropy_laplacian_trace_geometry_degeneracy.py"
OUT = ROOT / "generated" / "p2774_s1724_entropy_laplacian_trace_geometry_degeneracy.json"
MD = ROOT / "generated" / "p2774_s1724_entropy_laplacian_trace_geometry_degeneracy.md"


class P2774EntropyLaplacianTraceGeometryDegeneracyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2774_ENTROPY_LAPLACIAN_TRACE_GEOMETRY_DEGENERACY_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2773"], "P2773_SHANNON_ENTROPY_GEOMETRY_FORCING_OBSTRUCTION_NO_CLOSURE")
        self.assertIn("Laplacian trace", self.payload["audited_question"])

    def test_same_entropy_laplacian_trace_but_different_metric(self):
        witness = self.payload["entropy_laplacian_degeneracy_witness"]
        self.assertTrue(witness["entropy_matches_alpha_geo"])
        self.assertTrue(witness["same_degree_sequence"])
        self.assertTrue(witness["same_laplacian_trace"])
        self.assertEqual(witness["geometry_row_count"], 2)
        self.assertEqual(witness["distinct_distance_histogram_count"], 2)
        self.assertFalse(witness["entropy_plus_laplacian_trace_forces_unique_geometry"])
        rows = {row["geometry"]: row for row in witness["geometry_rows"]}
        self.assertEqual(set(rows), {"torus_4x4", "circulant_pm1_pm2"})
        self.assertEqual(rows["torus_4x4"]["degree_sequence"], rows["circulant_pm1_pm2"]["degree_sequence"])
        self.assertEqual(rows["torus_4x4"]["laplacian_trace"], rows["circulant_pm1_pm2"]["laplacian_trace"])
        self.assertNotEqual(rows["torus_4x4"]["distance_histogram"], rows["circulant_pm1_pm2"]["distance_histogram"])

    def test_acceptance_blocks_entropy_laplacian_geometry_forcing(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["facts"]["extra_laplacian_trace_principle_supplied"])
        self.assertTrue(acceptance["facts"]["same_degree_sequence_and_laplacian_trace"])
        self.assertFalse(acceptance["facts"]["entropy_plus_laplacian_trace_forces_unique_geometry"])
        self.assertFalse(acceptance["accepted_as_entropy_laplacian_geometry_forcing_theorem"])
        self.assertFalse(acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(acceptance["accepted_as_kernel_full_expression_theorem"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("entropy_plus_laplacian_trace_forces_unique_geometry", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("full Laplacian-spectrum", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2774", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2774/S1724", MD.read_text(encoding="utf-8"))
        self.assertIn("P2774/S1724", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2774/S1724", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2774", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
