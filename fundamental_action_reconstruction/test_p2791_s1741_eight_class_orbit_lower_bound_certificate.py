import json
import math
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
GEN = ROOT / "generated"
JSON_PATH = GEN / "p2791_s1741_eight_class_orbit_lower_bound_certificate.json"
MD_PATH = GEN / "p2791_s1741_eight_class_orbit_lower_bound_certificate.md"


class P2791EightClassOrbitLowerBoundCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON_PATH.exists():
            import subprocess
            subprocess.run(["python", str(ROOT / "p2791_s1741_eight_class_orbit_lower_bound_certificate.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.witness = cls.payload["orbit_lower_bound_witness"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2791_EIGHT_CLASS_ORBIT_LOWER_BOUND_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2786"], "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2790"], "P2790_EIGHTH_16NODE_WITNESS_NO_EXHAUSTION_CERTIFICATE_NO_CLOSURE")

    def test_pairwise_nonisomorphism_and_lower_bound(self):
        self.assertEqual(self.witness["representative_count"], 8)
        self.assertEqual(self.witness["pair_count"], 28)
        self.assertTrue(self.witness["all_28_pairs_nonisomorphic"])
        self.assertEqual(self.witness["certified_disjoint_labeled_orbit_lower_bound"], 13463256807000)
        self.assertEqual(self.witness["largest_single_orbit_label"], "p2790_eighth_witness")
        self.assertTrue(all(not row["isomorphic"] for row in self.witness["pairwise_isomorphism_rows"]))

    def test_orbit_stabilizer_rows_are_exact(self):
        expected_stabilizers = {
            "circulant_pm1_pm2": 32,
            "circulant_pm1_pm3": 32,
            "circulant_pm1_pm4": 32,
            "circulant_pm1_pm6": 32,
            "circulant_pm1_pm7": 4096,
            "torus_4x4": 384,
            "two_c8_layers_cross_pm0_pm4": 64,
            "p2790_eighth_witness": 2,
        }
        orbit_sum = 0
        for row in self.witness["representative_rows"]:
            self.assertEqual(row["edge_count"], 32)
            self.assertEqual(row["automorphism_group_size"], expected_stabilizers[row["label"]])
            self.assertEqual(math.factorial(16) // row["automorphism_group_size"], row["orbit_size_by_orbit_stabilizer"])
            orbit_sum += row["orbit_size_by_orbit_stabilizer"]
        self.assertEqual(orbit_sum, self.witness["certified_disjoint_labeled_orbit_lower_bound"])

    def test_acceptance_blocks_closure(self):
        self.assertTrue(self.acceptance["accepted_as_eight_class_orbit_lower_bound_certificate"])
        self.assertFalse(self.acceptance["accepted_as_full_16node_canonical_generator_certificate"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("P2697-P2791", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_updated(self):
        self.assertIn("P2791/S1741", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2791/S1741", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2791/S1741", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2791", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
