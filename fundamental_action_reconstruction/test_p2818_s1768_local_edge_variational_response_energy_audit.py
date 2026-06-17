import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2818_s1768_local_edge_variational_response_energy_audit.py"
JSON_PATH = ROOT / "generated" / "p2818_s1768_local_edge_variational_response_energy_audit.json"
MD_PATH = ROOT / "generated" / "p2818_s1768_local_edge_variational_response_energy_audit.md"


class P2818LocalEdgeVariationalResponseEnergyAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2818_LOCAL_EDGE_VARIATIONAL_RESPONSE_ENERGY_OBSTRUCTED_NO_COUPLING_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2817"], "P2817_STRUCTURAL_OBSERVABLE_SOURCE_CANDIDATE_OBSTRUCTED_NO_COUPLING_NO_CLOSURE")

    def test_local_edge_response_counts(self):
        audit = self.payload["e_local_source_candidate_audit"]
        self.assertEqual(audit["decoded_graph_count"], 16828)
        self.assertEqual(audit["e_local_class_count"], 344)
        self.assertEqual(audit["e_local_collision_class_count"], 272)
        self.assertEqual(audit["e_local_collision_graph_count"], 16756)
        self.assertEqual(audit["e_local_max_class_size"], 788)
        self.assertEqual(audit["remaining_defect_canonical_minus_e_local"], 16484)

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        self.assertTrue(matrix["facts"]["p2817_structural_observable_rejected"])
        self.assertTrue(matrix["facts"]["exactly_one_local_edge_response_energy_tested"])
        self.assertTrue(matrix["facts"]["energy_has_explicit_dimensionless_normalization_candidate"])
        self.assertFalse(matrix["facts"]["energy_separates_full_carrier"])
        self.assertFalse(matrix["facts"]["strict_graph_source_law_exported"])
        self.assertFalse(matrix["facts"]["typed_graph_to_K_or_Ltotal_coupling_exported"])
        self.assertTrue(matrix["accepted_as_local_edge_response_energy_audit"])
        self.assertFalse(matrix["accepted_as_complete_carrier_separator"])
        self.assertFalse(matrix["accepted_as_source_law_coupling"])
        self.assertTrue(matrix["accepted_as_bounded_candidate_rejection"])

    def test_documents_updated(self):
        self.assertIn("P2818/S1768", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2818/S1768", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2818/S1768", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2818", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
