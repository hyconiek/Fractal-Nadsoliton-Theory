import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2833_s1783_edge_toggle_response_polynomial_source_audit.py"
JSON_PATH = ROOT / "generated" / "p2833_s1783_edge_toggle_response_polynomial_source_audit.json"
MD_PATH = ROOT / "generated" / "p2833_s1783_edge_toggle_response_polynomial_source_audit.md"


class P2833EdgeToggleResponsePolynomialSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(
            self.payload["status"],
            "P2833_EDGE_TOGGLE_RESPONSE_POLYNOMIAL_RESIDUAL_NO_GO_NO_COUPLING_NO_CLOSURE",
        )
        self.assertEqual(
            self.payload["input_statuses"]["P2818"],
            "P2818_LOCAL_EDGE_VARIATIONAL_RESPONSE_ENERGY_OBSTRUCTED_NO_COUPLING_NO_CLOSURE",
        )
        self.assertEqual(
            self.payload["input_statuses"]["P2832"],
            "P2832_COMMON_NEIGHBOR_CYCLE_SOURCE_FORMULA_NO_GO_NO_COUPLING_NO_CLOSURE",
        )

    def test_edge_toggle_counts(self):
        audit = self.payload["edge_toggle_response_polynomial_source_audit"]
        self.assertEqual(audit["decoded_graph_count"], 16828)
        self.assertTrue(audit["higher_order_non_label_graph_formula_exported_for_candidate"])
        self.assertFalse(audit["digest_label_hash_or_rank_used"])
        self.assertTrue(audit["edge_toggle_variational_response_coefficients_computed"])
        self.assertTrue(audit["dimensionless_normalization_available_for_candidate"])
        self.assertEqual(audit["candidate_class_count"], 16757)
        self.assertEqual(audit["candidate_collision_class_count"], 67)
        self.assertEqual(audit["candidate_collision_graph_count"], 138)
        self.assertEqual(audit["candidate_max_class_size"], 3)
        self.assertEqual(audit["candidate_defect_after_formula"], 71)
        self.assertEqual(audit["candidate_class_size_histogram"]["1"], 16690)
        self.assertEqual(audit["candidate_class_size_histogram"]["2"], 63)
        self.assertEqual(audit["candidate_class_size_histogram"]["3"], 4)

    def test_largest_collision_samples(self):
        largest = self.payload["edge_toggle_response_polynomial_source_audit"]["largest_collision_samples"]
        self.assertEqual(largest[0]["class_size"], 3)
        self.assertGreaterEqual(len(largest[0]["sample_graph_indices_0_based"]), 3)
        self.assertIsInstance(largest[0]["response_polynomial"], list)

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        facts = matrix["facts"]
        self.assertTrue(facts["p2818_local_response_rejected"])
        self.assertTrue(facts["p2832_low_order_formula_rejected"])
        self.assertTrue(facts["exactly_one_higher_order_non_label_formula_tested"])
        self.assertFalse(facts["digest_label_hash_or_rank_used"])
        self.assertTrue(facts["edge_toggle_variational_response_coefficients_computed"])
        self.assertTrue(facts["dimensionless_normalization_available_for_candidate"])
        self.assertFalse(facts["candidate_separates_full_carrier"])
        self.assertFalse(facts["proved_variational_derivative_into_K_or_Ltotal_exported"])
        self.assertFalse(facts["typed_graph_to_K_or_Ltotal_coupling_theorem_exported"])
        self.assertFalse(matrix["accepted_as_source_law_coupling"])
        self.assertTrue(matrix["accepted_as_bounded_edge_toggle_witness_with_residual_no_go"])

    def test_documents_updated(self):
        self.assertIn("P2833/S1783", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2833/S1783", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2833/S1783", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2833", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
