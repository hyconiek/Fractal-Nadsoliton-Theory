import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2819_s1769_p2815_digest_blocker_cut_response_audit.py"
JSON_PATH = ROOT / "generated" / "p2819_s1769_p2815_digest_blocker_cut_response_audit.json"
MD_PATH = ROOT / "generated" / "p2819_s1769_p2815_digest_blocker_cut_response_audit.md"


class P2819P2815DigestBlockerCutResponseAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2819_P2815_DIGEST_BLOCKER_CUT_SAMPLE_WITNESS_NO_FULL_COUPLING_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2818"], "P2818_LOCAL_EDGE_VARIATIONAL_RESPONSE_ENERGY_OBSTRUCTED_NO_COUPLING_NO_CLOSURE")

    def test_sample_response_counts(self):
        audit = self.payload["p2815_digest_blocker_cut_audit"]
        self.assertEqual(audit["decoded_graph_count"], 16828)
        self.assertEqual(audit["p2818_blocker_class_size"], 788)
        self.assertEqual(audit["computed_graph_count"], 24)
        self.assertEqual(audit["sample_fraction_of_blocker_cut"], [24, 788])
        self.assertEqual(audit["toggle_digest_refined_class_count"], 24)
        self.assertEqual(audit["toggle_digest_collision_class_count"], 0)
        self.assertEqual(audit["toggle_digest_collision_graph_count"], 0)
        self.assertEqual(audit["toggle_digest_max_class_size"], 1)
        self.assertEqual(audit["sample_defect_after_toggle_digest"], 0)

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        self.assertTrue(matrix["facts"]["p2818_local_response_rejected"])
        self.assertTrue(matrix["facts"]["exactly_one_blocker_cut_tested"])
        self.assertTrue(matrix["facts"]["uses_full_p2815_spectral_local_motif_response_digest"])
        self.assertTrue(matrix["facts"]["blocker_cut_is_fully_separated"])
        self.assertFalse(matrix["facts"]["all_16828_graphs_audited"])
        self.assertFalse(matrix["facts"]["entire_blocker_cut_audited"])
        self.assertFalse(matrix["facts"]["strict_graph_source_law_exported"])
        self.assertFalse(matrix["facts"]["typed_graph_to_K_or_Ltotal_coupling_exported"])
        self.assertTrue(matrix["accepted_as_sample_response_witness"])
        self.assertFalse(matrix["accepted_as_full_carrier_source_law_coupling"])
        self.assertTrue(matrix["accepted_as_bounded_no_closure_audit"])

    def test_documents_updated(self):
        self.assertIn("P2819/S1769", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2819/S1769", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2819/S1769", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2819", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
