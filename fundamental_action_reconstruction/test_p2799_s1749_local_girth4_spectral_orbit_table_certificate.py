import json
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
GEN = ROOT / "generated"
JSON_PATH = GEN / "p2799_s1749_local_girth4_spectral_orbit_table_certificate.json"
MD_PATH = GEN / "p2799_s1749_local_girth4_spectral_orbit_table_certificate.md"


class P2799LocalGirth4SpectralOrbitTableCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON_PATH.exists():
            import subprocess
            subprocess.run(["python", str(ROOT / "p2799_s1749_local_girth4_spectral_orbit_table_certificate.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.witness = cls.payload["local_girth4_spectral_orbit_witness"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2799_LOCAL_GIRTH4_SPECTRAL_ORBIT_TABLE_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2786"], "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2790"], "P2790_EIGHTH_16NODE_WITNESS_NO_EXHAUSTION_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2798"], "P2798_EXTERNAL_GIRTH4_SUBTARGET_LOCAL_GIRTH_FILTER_CERTIFICATE_NO_CLOSURE")

    def test_table_summary(self):
        self.assertEqual(self.witness["input_girth4_label_count"], 6)
        self.assertEqual(self.witness["table_row_count"], 6)
        self.assertTrue(self.witness["all_rows_16node_4regular"])
        self.assertTrue(self.witness["all_complements_16node_11regular"])
        self.assertTrue(self.witness["all_automorphism_counts_match_complements"])
        self.assertEqual(self.witness["distinct_adjacency_charpoly_count"], 6)
        self.assertEqual(self.witness["total_labeled_orbit_size_for_local_girth4_rows"], 2348024679000)

    def test_each_row_has_exact_polynomial_and_orbit_data(self):
        labels = [row["label"] for row in self.witness["table_rows"]]
        self.assertEqual(labels, [
            "circulant_pm1_pm3",
            "circulant_pm1_pm4",
            "circulant_pm1_pm6",
            "circulant_pm1_pm7",
            "torus_4x4",
            "two_c8_layers_cross_pm0_pm4",
        ])
        for row in self.witness["table_rows"]:
            self.assertEqual(row["edge_count"], 32)
            self.assertEqual(row["degree_sequence"], [4] * 16)
            self.assertEqual(row["complement_edge_count"], 88)
            self.assertEqual(row["complement_degree_sequence"], [11] * 16)
            self.assertTrue(row["automorphism_matches_complement"])
            self.assertEqual(len(row["adjacency_charpoly_coefficients"]), 17)
            self.assertEqual(len(row["laplacian_charpoly_coefficients"]), 17)
            self.assertEqual(len(row["signless_laplacian_charpoly_coefficients"]), 17)
            self.assertEqual(len(row["complement_adjacency_charpoly_coefficients"]), 17)
            self.assertEqual(len(row["complement_laplacian_charpoly_coefficients"]), 17)
            self.assertEqual(len(row["complement_signless_laplacian_charpoly_coefficients"]), 17)

    def test_acceptance_blocks_closure(self):
        self.assertTrue(self.acceptance["accepted_as_local_girth4_spectral_orbit_table_certificate"])
        self.assertFalse(self.acceptance["accepted_as_full_16node_canonical_generator_certificate"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("P2697-P2799", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_updated(self):
        self.assertIn("P2799/S1749", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2799/S1749", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2799/S1749", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2799", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
