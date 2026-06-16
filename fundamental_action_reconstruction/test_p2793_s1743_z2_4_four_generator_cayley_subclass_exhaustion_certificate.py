import json
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
GEN = ROOT / "generated"
JSON_PATH = GEN / "p2793_s1743_z2_4_four_generator_cayley_subclass_exhaustion_certificate.json"
MD_PATH = GEN / "p2793_s1743_z2_4_four_generator_cayley_subclass_exhaustion_certificate.md"


class P2793Z2FourGeneratorCayleySubclassExhaustionCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON_PATH.exists():
            import subprocess
            subprocess.run(["python", str(ROOT / "p2793_s1743_z2_4_four_generator_cayley_subclass_exhaustion_certificate.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.witness = cls.payload["z2_4_cayley_subclass_witness"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2793_Z2_4_FOUR_GENERATOR_CAYLEY_SUBCLASS_EXHAUSTION_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2786"], "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2791"], "P2791_EIGHT_CLASS_ORBIT_LOWER_BOUND_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2792"], "P2792_C16_TWO_JUMP_CIRCULANT_SUBCLASS_EXHAUSTION_CERTIFICATE_NO_CLOSURE")

    def test_complete_z2_4_cayley_subclass_counts(self):
        self.assertEqual(self.witness["raw_generator_set_count"], 1365)
        self.assertEqual(self.witness["connected_generator_set_count"], 840)
        self.assertEqual(self.witness["disconnected_generator_set_count"], 525)
        self.assertEqual(self.witness["rank_histogram"], {"1": 0, "2": 0, "3": 525, "4": 840})
        self.assertEqual(self.witness["connected_isomorphism_class_count"], 1)
        self.assertEqual(self.witness["pair_count_after_quotient"], 0)
        self.assertTrue(self.witness["all_connected_classes_pairwise_nonisomorphic"])
        self.assertTrue(self.witness["all_connected_classes_represented_in_p2786_local_set"])
        self.assertTrue(self.witness["all_connected_classes_represented_in_p2791_eight_class_set"])

    def test_single_class_matches_known_torus(self):
        self.assertEqual(len(self.witness["class_rows"]), 1)
        row = self.witness["class_rows"][0]
        self.assertEqual(row["representative_generators"], [1, 2, 4, 8])
        self.assertEqual(row["generator_set_member_count"], 840)
        self.assertEqual(row["matching_p2786_local_representatives"], ["torus_4x4"])
        self.assertEqual(row["matching_p2791_eight_class_representatives"], ["torus_4x4"])
        self.assertEqual(row["automorphism_group_size"], 384)
        self.assertEqual(row["orbit_size_by_orbit_stabilizer"], 54486432000)
        self.assertEqual(len(row["adjacency_charpoly_coefficients"]), 17)

    def test_acceptance_blocks_closure(self):
        self.assertTrue(self.acceptance["accepted_as_complete_z2_4_four_generator_cayley_subclass_certificate"])
        self.assertFalse(self.acceptance["accepted_as_full_16node_canonical_generator_certificate"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("P2697-P2793", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_updated(self):
        self.assertIn("P2793/S1743", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2793/S1743", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2793/S1743", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2793", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
