import json
import math
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
GEN = ROOT / "generated"
JSON_PATH = GEN / "p2789_s1739_orbit_stabilizer_exact_quotient_certificate.json"
MD_PATH = GEN / "p2789_s1739_orbit_stabilizer_exact_quotient_certificate.md"


class P2789OrbitStabilizerExactQuotientCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON_PATH.exists():
            import subprocess
            subprocess.run(["python", str(ROOT / "p2789_s1739_orbit_stabilizer_exact_quotient_certificate.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.witness = cls.payload["orbit_stabilizer_witness"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2789_ORBIT_STABILIZER_EXACT_QUOTIENT_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2786"], "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2787"], "P2787_SMALL_CANONICAL_GENERATOR_PIPELINE_AUDIT_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2788"], "P2788_COMPLEMENT_DUALITY_EXACT_SPECTRAL_CERTIFICATE_NO_CLOSURE")

    def test_complete_small_orbit_stabilizer_arithmetic(self):
        self.assertEqual(self.witness["small_complete_8node_row_count"], 6)
        self.assertEqual(self.witness["small_orbit_size_sum"], 19355)
        self.assertEqual(self.witness["stored_small_connected_labeled_candidate_count"], 19355)
        self.assertTrue(self.witness["all_small_orbit_stabilizer_counts_match_stored_members"])
        self.assertTrue(self.witness["small_orbit_sum_matches_stored_connected_labeled_total"])
        for row in self.witness["small_8node_rows"]:
            self.assertEqual(math.factorial(8) // row["automorphism_group_size"], row["stored_labeled_member_count"])
            self.assertTrue(row["orbit_stabilizer_matches_stored_member_count"])

    def test_local_16node_stabilizers_are_exact_positive_divisors(self):
        expected_stabilizers = {
            "circulant_pm1_pm2": 32,
            "circulant_pm1_pm3": 32,
            "circulant_pm1_pm4": 32,
            "circulant_pm1_pm6": 32,
            "circulant_pm1_pm7": 4096,
            "torus_4x4": 384,
            "two_c8_layers_cross_pm0_pm4": 64,
        }
        self.assertEqual(self.witness["local_16node_row_count"], 7)
        self.assertTrue(self.witness["all_local_16node_stabilizers_positive"])
        for row in self.witness["local_16node_rows"]:
            self.assertEqual(row["n"], 16)
            self.assertEqual(row["edge_count"], 32)
            self.assertEqual(row["automorphism_group_size"], expected_stabilizers[row["representative"]])
            self.assertEqual(math.factorial(16) % row["automorphism_group_size"], 0)
            self.assertEqual(math.factorial(16) // row["automorphism_group_size"], row["orbit_size_by_orbit_stabilizer"])

    def test_acceptance_blocks_closure(self):
        self.assertTrue(self.acceptance["accepted_as_exact_orbit_stabilizer_quotient_certificate"])
        self.assertFalse(self.acceptance["accepted_as_full_16node_canonical_generator_certificate"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("P2697-P2789", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_updated(self):
        self.assertIn("P2789/S1739", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2789/S1739", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2789/S1739", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2789", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
