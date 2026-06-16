import json
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
GEN = ROOT / "generated"
JSON_PATH = GEN / "p2792_s1742_c16_two_jump_circulant_subclass_exhaustion_certificate.json"
MD_PATH = GEN / "p2792_s1742_c16_two_jump_circulant_subclass_exhaustion_certificate.md"


class P2792C16TwoJumpCirculantSubclassExhaustionCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON_PATH.exists():
            import subprocess
            subprocess.run(["python", str(ROOT / "p2792_s1742_c16_two_jump_circulant_subclass_exhaustion_certificate.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.witness = cls.payload["circulant_subclass_witness"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2792_C16_TWO_JUMP_CIRCULANT_SUBCLASS_EXHAUSTION_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2786"], "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2791"], "P2791_EIGHT_CLASS_ORBIT_LOWER_BOUND_CERTIFICATE_NO_CLOSURE")

    def test_complete_named_subclass_counts(self):
        self.assertEqual(self.witness["raw_jump_pair_count"], 21)
        self.assertEqual(self.witness["connected_jump_pair_count"], 18)
        self.assertEqual(self.witness["disconnected_jump_pair_count"], 3)
        self.assertEqual(self.witness["connected_isomorphism_class_count"], 5)
        self.assertEqual(self.witness["pair_count_after_quotient"], 10)
        self.assertTrue(self.witness["all_connected_classes_pairwise_nonisomorphic"])
        self.assertTrue(self.witness["all_connected_classes_represented_in_p2786_local_set"])
        self.assertTrue(all(not row["isomorphic"] for row in self.witness["pairwise_class_rows"]))

    def test_class_rows_have_expected_local_matches(self):
        expected_matches = [
            ["circulant_pm1_pm2"],
            ["circulant_pm1_pm3"],
            ["circulant_pm1_pm4"],
            ["circulant_pm1_pm6"],
            ["circulant_pm1_pm7"],
        ]
        self.assertEqual([row["matching_p2786_local_representatives"] for row in self.witness["class_rows"]], expected_matches)
        self.assertEqual([row["jump_pair_member_count"] for row in self.witness["class_rows"]], [4, 4, 4, 4, 2])
        for row in self.witness["class_rows"]:
            self.assertEqual(row["edge_count"], 32)
            self.assertEqual(len(row["adjacency_charpoly_coefficients"]), 17)
            self.assertTrue(row["represented_in_p2786_local_set"])

    def test_acceptance_blocks_closure(self):
        self.assertTrue(self.acceptance["accepted_as_complete_c16_two_jump_circulant_subclass_certificate"])
        self.assertFalse(self.acceptance["accepted_as_full_16node_canonical_generator_certificate"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("P2697-P2792", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_updated(self):
        self.assertIn("P2792/S1742", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2792/S1742", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2792/S1742", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2792", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
