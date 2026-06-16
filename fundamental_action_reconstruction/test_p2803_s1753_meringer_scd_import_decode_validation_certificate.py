import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
JSON_PATH = ROOT / "generated" / "p2803_s1753_meringer_scd_import_decode_validation_certificate.json"
MD_PATH = ROOT / "generated" / "p2803_s1753_meringer_scd_import_decode_validation_certificate.md"


class P2803MeringerScdImportDecodeValidationCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(ROOT / "p2803_s1753_meringer_scd_import_decode_validation_certificate.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.cert = cls.payload["meringer_scd_import_decode_validation"]
        cls.parse = cls.cert["parse_stats"]
        cls.metrics = cls.cert["graph_metrics"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_prior_obstruction_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2803_MERINGER_SCD_IMPORT_DECODE_VALIDATION_CERTIFICATE_IMPORT_VALIDATED_NO_QUOTIENT_NO_CLOSURE",
        )
        self.assertEqual(
            self.payload["input_statuses"]["P2802"],
            "P2802_GIRTH4_FETCH_OBSTRUCTION_TAXONOMY_CERTIFICATE_NO_IMPORT_NO_CLOSURE",
        )

    def test_artifact_hash_and_byte_size(self):
        self.assertTrue(self.cert["file_present"])
        self.assertTrue(self.cert["decoder_present"])
        self.assertEqual(self.cert["byte_size"], 150489)
        self.assertEqual(
            self.cert["sha256"],
            "160bf01bc0767652bb05c0466a9358628fd5e5053672695309a04e307fe25a88",
        )
        self.assertTrue(self.cert["sha256_matches_expected"])

    def test_scd_parse_count_and_integrity(self):
        self.assertEqual(self.parse["encoded_edge_code_length_bytes"], 32)
        self.assertEqual(self.parse["decoded_entry_count"], 16828)
        self.assertTrue(self.parse["parse_consumed_all_bytes"])
        self.assertEqual(self.parse["malformed_prefix_count"], 0)
        self.assertEqual(self.parse["truncated_tail_count"], 0)
        self.assertEqual(self.parse["code_length_mismatch_count"], 0)

    def test_decoded_graph_invariants(self):
        self.assertEqual(self.metrics["strict_4_regular_graph_count"], 16828)
        self.assertEqual(self.metrics["girth_at_least_4_no_triangle_graph_count"], 16828)
        self.assertEqual(self.metrics["connected_graph_count"], 16828)
        self.assertEqual(self.metrics["bad_4_regular_graph_count"], 0)
        self.assertEqual(self.metrics["triangle_graph_count"], 0)
        self.assertEqual(self.metrics["disconnected_graph_count"], 0)
        self.assertEqual(self.metrics["simple_loop_violation_count"], 0)
        self.assertEqual(self.metrics["duplicate_neighbor_violation_count"], 0)
        self.assertEqual(self.metrics["symmetric_adjacency_violation_count"], 0)
        self.assertEqual(self.metrics["edge_count_histogram"], {"32": 16828})
        self.assertEqual(self.metrics["unique_decoded_adjacency_hash_count"], 16828)

    def test_acceptance_boundaries(self):
        self.assertTrue(self.acceptance["accepted_as_meringer_scd_import_decode_validation"])
        self.assertFalse(self.acceptance["accepted_as_completed_girth4_quotient_audit"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("quotient_charpoly_complement_orbit_audit_completed", self.acceptance["missing_criteria"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["toe_closure_exported"])

    def test_written_documents_reference_guardrail(self):
        self.assertIn("P2803/S1753", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2803/S1753", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2803/S1753", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2803", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
