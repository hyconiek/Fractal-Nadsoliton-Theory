import json
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
GEN = ROOT / "generated"
JSON_PATH = GEN / "p2798_s1748_external_girth4_subtarget_local_girth_filter_certificate.json"
MD_PATH = GEN / "p2798_s1748_external_girth4_subtarget_local_girth_filter_certificate.md"


class P2798ExternalGirth4SubtargetLocalGirthFilterCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON_PATH.exists():
            import subprocess
            subprocess.run(["python", str(ROOT / "p2798_s1748_external_girth4_subtarget_local_girth_filter_certificate.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.witness = cls.payload["external_girth4_subtarget_witness"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2798_EXTERNAL_GIRTH4_SUBTARGET_LOCAL_GIRTH_FILTER_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2786"], "P2786_GRAPH6_PROVENANCE_TOOLCHAIN_GATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2790"], "P2790_EIGHTH_16NODE_WITNESS_NO_EXHAUSTION_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2791"], "P2791_EIGHT_CLASS_ORBIT_LOWER_BOUND_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2797"], "P2797_EXTERNAL_FULL_COUNT_GAP_CERTIFICATE_NO_CLOSURE")

    def test_external_girth4_gap_numbers(self):
        self.assertEqual(self.witness["external_connected_16node_4regular_girth4_class_count"], 16828)
        self.assertEqual(self.witness["local_representative_count"], 8)
        self.assertEqual(self.witness["local_girth_at_least_4_count"], 6)
        self.assertEqual(self.witness["local_girth_at_least_4_labels"], [
            "circulant_pm1_pm3",
            "circulant_pm1_pm4",
            "circulant_pm1_pm6",
            "circulant_pm1_pm7",
            "torus_4x4",
            "two_c8_layers_cross_pm0_pm4",
        ])
        self.assertEqual(self.witness["local_triangle_containing_labels"], ["circulant_pm1_pm2", "p2790_eighth_witness"])
        self.assertEqual(self.witness["external_girth4_gap_after_current_witnesses"], 16822)
        self.assertEqual(self.witness["local_girth4_coverage_fraction_exact"], "3/8414")

    def test_local_girth_rows_are_exactly_eight(self):
        rows = self.witness["local_representative_girth_rows"]
        self.assertEqual(len(rows), 8)
        self.assertEqual([row["girth"] for row in rows], [3, 4, 4, 4, 4, 4, 4, 3])
        self.assertTrue(all(row["edge_count"] == 32 for row in rows))

    def test_external_source_is_scoped(self):
        source = self.witness["external_source"]
        self.assertEqual(source["url"], "https://www.mathe2.uni-bayreuth.de/markus/reggraphs.html")
        self.assertEqual(source["linked_detail_url"], "https://www.mathe2.uni-bayreuth.de/markus/REGGRAPHS/16_4_4.html")
        self.assertEqual(source["table"], "Connected regular graphs with girth at least 4")
        self.assertEqual(source["reported_count"], 16828)

    def test_acceptance_blocks_closure(self):
        self.assertTrue(self.acceptance["accepted_as_external_girth4_subtarget_local_girth_filter_certificate"])
        self.assertFalse(self.acceptance["accepted_as_full_16node_canonical_generator_certificate"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("P2697-P2798", self.payload["decision"]["next_honest_step"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_updated(self):
        self.assertIn("P2798/S1748", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2798/S1748", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2798/S1748", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2798", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
