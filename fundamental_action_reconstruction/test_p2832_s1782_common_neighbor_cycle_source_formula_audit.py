import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2832_s1782_common_neighbor_cycle_source_formula_audit.py"
JSON_PATH = ROOT / "generated" / "p2832_s1782_common_neighbor_cycle_source_formula_audit.json"
MD_PATH = ROOT / "generated" / "p2832_s1782_common_neighbor_cycle_source_formula_audit.md"


class P2832CommonNeighborCycleSourceFormulaAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2832_COMMON_NEIGHBOR_CYCLE_SOURCE_FORMULA_NO_GO_NO_COUPLING_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2818"], "P2818_LOCAL_EDGE_VARIATIONAL_RESPONSE_ENERGY_OBSTRUCTED_NO_COUPLING_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2831"], "P2831_P2815_DIGEST_SOURCE_LAW_COUPLING_OBLIGATION_NO_GO_NO_CLOSURE")

    def test_formula_counts(self):
        audit = self.payload["common_neighbor_cycle_source_formula_audit"]
        self.assertEqual(audit["decoded_graph_count"], 16828)
        self.assertTrue(audit["non_label_graph_formula_exported_for_candidate"])
        self.assertTrue(audit["dimensionless_normalization_available_for_candidate"])
        self.assertEqual(audit["candidate_class_count"], 344)
        self.assertEqual(audit["candidate_collision_class_count"], 272)
        self.assertEqual(audit["candidate_collision_graph_count"], 16756)
        self.assertEqual(audit["candidate_max_class_size"], 788)
        self.assertEqual(audit["candidate_defect_after_formula"], 16484)
        self.assertEqual(audit["candidate_class_size_histogram"]["1"], 72)
        self.assertEqual(audit["candidate_class_size_histogram"]["788"], 1)

    def test_largest_collision_samples(self):
        audit = self.payload["common_neighbor_cycle_source_formula_audit"]
        largest = audit["largest_collision_samples"]
        self.assertEqual(largest[0]["class_size"], 788)
        self.assertEqual(len(largest[0]["sample_graph_indices_0_based"]), 12)
        self.assertEqual(largest[0]["profile"][0], 0)
        self.assertIsInstance(largest[0]["profile"][2], list)

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        self.assertTrue(matrix["facts"]["p2818_local_response_rejected"])
        self.assertTrue(matrix["facts"]["p2831_digest_label_lane_rejected"])
        self.assertTrue(matrix["facts"]["exactly_one_non_label_formula_tested"])
        self.assertTrue(matrix["facts"]["non_label_graph_formula_exported_for_candidate"])
        self.assertTrue(matrix["facts"]["dimensionless_normalization_available_for_candidate"])
        self.assertFalse(matrix["facts"]["candidate_separates_full_carrier"])
        self.assertFalse(matrix["facts"]["variational_derivative_exported"])
        self.assertFalse(matrix["facts"]["typed_graph_to_K_or_Ltotal_coupling_theorem_exported"])
        self.assertFalse(matrix["accepted_as_source_law_coupling"])
        self.assertTrue(matrix["accepted_as_bounded_formula_no_go"])

    def test_documents_updated(self):
        self.assertIn("P2832/S1782", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2832/S1782", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2832/S1782", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2832", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
