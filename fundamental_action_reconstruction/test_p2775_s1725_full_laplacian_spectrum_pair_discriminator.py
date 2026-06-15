import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2775_s1725_full_laplacian_spectrum_pair_discriminator.py"
OUT = ROOT / "generated" / "p2775_s1725_full_laplacian_spectrum_pair_discriminator.json"
MD = ROOT / "generated" / "p2775_s1725_full_laplacian_spectrum_pair_discriminator.md"


class P2775FullLaplacianSpectrumPairDiscriminatorTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2775_FULL_LAPLACIAN_SPECTRUM_PAIR_DISCRIMINATOR_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2774"], "P2774_ENTROPY_LAPLACIAN_TRACE_GEOMETRY_DEGENERACY_NO_CLOSURE")
        self.assertIn("full Laplacian spectrum", self.payload["audited_question"])

    def test_full_spectrum_breaks_exact_p2774_pair_degeneracy(self):
        witness = self.payload["full_laplacian_spectrum_witness"]
        self.assertTrue(witness["entropy_matches_alpha_geo"])
        self.assertTrue(witness["same_laplacian_trace"])
        self.assertEqual(witness["geometry_row_count"], 2)
        self.assertEqual(witness["distinct_full_laplacian_spectrum_count"], 2)
        self.assertTrue(witness["full_spectrum_breaks_p2774_pair_degeneracy"])
        rows = {row["geometry"]: row for row in witness["geometry_rows"]}
        self.assertEqual(set(rows), {"torus_4x4", "circulant_pm1_pm2"})
        self.assertEqual(rows["torus_4x4"]["laplacian_trace"], rows["circulant_pm1_pm2"]["laplacian_trace"])
        self.assertNotEqual(rows["torus_4x4"]["laplacian_spectrum_rounded"], rows["circulant_pm1_pm2"]["laplacian_spectrum_rounded"])
        self.assertAlmostEqual(rows["torus_4x4"]["spectral_sum"], rows["circulant_pm1_pm2"]["spectral_sum"], places=8)
        self.assertEqual(rows["torus_4x4"]["spectral_zero_multiplicity"], 1)
        self.assertEqual(rows["circulant_pm1_pm2"]["spectral_zero_multiplicity"], 1)

    def test_acceptance_is_positive_pair_local_but_blocks_closure(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["accepted_as_pair_local_degeneracy_breaking_witness"])
        self.assertFalse(acceptance["facts"]["sourced_nadsoliton_spectral_law_exported"])
        self.assertFalse(acceptance["facts"]["global_graph_class_uniqueness_audited"])
        self.assertFalse(acceptance["facts"]["kernel_or_ltotal_variational_coupling_exported"])
        self.assertFalse(acceptance["accepted_as_full_spectrum_nadsoliton_geometry_theorem"])
        self.assertFalse(acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("sourced_nadsoliton_spectral_law_exported", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("cospectral-degeneracy audit", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2775", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2775/S1725", MD.read_text(encoding="utf-8"))
        self.assertIn("P2775/S1725", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2775/S1725", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2775", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
