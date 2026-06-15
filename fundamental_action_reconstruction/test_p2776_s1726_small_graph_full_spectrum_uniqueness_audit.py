import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2776_s1726_small_graph_full_spectrum_uniqueness_audit.py"
OUT = ROOT / "generated" / "p2776_s1726_small_graph_full_spectrum_uniqueness_audit.json"
MD = ROOT / "generated" / "p2776_s1726_small_graph_full_spectrum_uniqueness_audit.md"


class P2776SmallGraphFullSpectrumUniquenessAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2776_SMALL_GRAPH_FULL_SPECTRUM_UNIQUENESS_AUDIT_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2775"], "P2775_FULL_LAPLACIAN_SPECTRUM_PAIR_DISCRIMINATOR_NO_CLOSURE")
        self.assertIn("full Laplacian spectrum", self.payload["audited_question"])

    def test_exhaustive_small_class_has_no_cospectral_collisions(self):
        witness = self.payload["small_graph_full_spectrum_uniqueness_witness"]
        self.assertEqual(witness["total_connected_unlabeled_classes"], 27)
        self.assertEqual(witness["total_cospectral_nonisomorphic_collision_count"], 0)
        self.assertTrue(witness["full_spectrum_injective_on_this_finite_class"])
        orders = {row["vertex_count"]: row for row in witness["orders"]}
        self.assertEqual(set(orders), {4, 5})
        self.assertEqual(orders[4]["connected_unlabeled_classes"], 6)
        self.assertEqual(orders[5]["connected_unlabeled_classes"], 21)
        self.assertEqual(orders[4]["distinct_laplacian_spectra"], 6)
        self.assertEqual(orders[5]["distinct_laplacian_spectra"], 21)
        self.assertTrue(orders[4]["full_spectrum_injective_on_unlabeled_connected_classes"])
        self.assertTrue(orders[5]["full_spectrum_injective_on_unlabeled_connected_classes"])

    def test_acceptance_is_small_class_positive_but_blocks_closure(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["accepted_as_small_class_spectral_uniqueness_witness"])
        self.assertFalse(acceptance["facts"]["strict_nadsoliton_spectral_source_law_exported"])
        self.assertFalse(acceptance["facts"]["sixteen_point_or_strict_graph_class_audited"])
        self.assertFalse(acceptance["facts"]["kernel_or_ltotal_variational_coupling_exported"])
        self.assertFalse(acceptance["accepted_as_global_full_spectrum_geometry_theorem"])
        self.assertFalse(acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("sixteen_point_or_strict_graph_class_audited", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("16-point/regular/strict graph class", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2776", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2776/S1726", MD.read_text(encoding="utf-8"))
        self.assertIn("P2776/S1726", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2776/S1726", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2776", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
