import json
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
GEN = ROOT / "generated"
JSON_PATH = GEN / "p2788_s1738_complement_duality_exact_spectral_certificate.json"
MD_PATH = GEN / "p2788_s1738_complement_duality_exact_spectral_certificate.md"


class P2788ComplementDualityExactSpectralCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON_PATH.exists():
            import subprocess
            subprocess.run(["python", str(ROOT / "p2788_s1738_complement_duality_exact_spectral_certificate.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.witness = cls.payload["complement_duality_witness"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2788_COMPLEMENT_DUALITY_EXACT_SPECTRAL_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2786"], "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2787"], "P2787_SMALL_CANONICAL_GENERATOR_PIPELINE_AUDIT_NO_CLOSURE")

    def test_exact_complement_duality_counts(self):
        self.assertEqual(self.witness["small_complete_8node_row_count"], 6)
        self.assertEqual(self.witness["local_16node_row_count"], 7)
        self.assertEqual(self.witness["total_rows_checked"], 13)
        self.assertTrue(self.witness["all_complements_regular_with_expected_degree"])
        self.assertTrue(self.witness["all_exact_adjacency_complement_identities_hold"])
        self.assertTrue(self.witness["all_exact_laplacian_complement_identities_hold"])
        self.assertTrue(self.witness["all_exact_signless_complement_identities_hold"])

    def test_complement_degrees_are_expected(self):
        for row in self.witness["small_complete_8node_rows"]:
            self.assertEqual(row["n"], 8)
            self.assertEqual(row["degree"], 4)
            self.assertEqual(row["complement_degree"], 3)
            self.assertTrue(row["all_exact_complement_identities_hold"])
        for row in self.witness["local_16node_rows"]:
            self.assertEqual(row["n"], 16)
            self.assertEqual(row["degree"], 4)
            self.assertEqual(row["complement_degree"], 11)
            self.assertTrue(row["all_exact_complement_identities_hold"])

    def test_acceptance_blocks_closure(self):
        self.assertTrue(self.acceptance["accepted_as_exact_complement_duality_certificate"])
        self.assertFalse(self.acceptance["accepted_as_full_16node_canonical_generator_certificate"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("P2697-P2788", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_updated(self):
        self.assertIn("P2788/S1738", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2788/S1738", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2788/S1738", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2788", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
