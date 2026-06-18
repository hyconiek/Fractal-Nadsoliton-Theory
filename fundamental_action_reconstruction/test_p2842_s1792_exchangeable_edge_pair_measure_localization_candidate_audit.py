import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2842_s1792_exchangeable_edge_pair_measure_localization_candidate_audit.py"
JSON_PATH = ROOT / "generated" / "p2842_s1792_exchangeable_edge_pair_measure_localization_candidate_audit.json"
MD_PATH = ROOT / "generated" / "p2842_s1792_exchangeable_edge_pair_measure_localization_candidate_audit.md"


class P2842ExchangeableEdgePairMeasureLocalizationCandidateAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status(self):
        self.assertEqual(
            self.payload["status"],
            "P2842_EXCHANGEABLE_EDGE_PAIR_MEASURE_LOCALIZATION_CANDIDATE_NO_GO_NO_COUPLING_NO_CLOSURE",
        )

    def test_candidate_measure(self):
        audit = self.payload["exchangeable_edge_pair_measure_localization_candidate_audit"]
        candidate = audit["candidate_measure"]
        self.assertEqual(candidate["candidate_object"], "mu_edge^exch(G)")
        self.assertEqual(candidate["vertex_count"], 16)
        self.assertEqual(candidate["regular_degree"], 4)
        self.assertEqual(candidate["unordered_pair_count"], 120)
        self.assertEqual(candidate["edge_count_per_graph"], 32)
        self.assertEqual(candidate["edge_density"]["numerator"], 4)
        self.assertEqual(candidate["edge_density"]["denominator"], 15)
        self.assertTrue(candidate["label_gauge_invariant_by_construction"])
        self.assertTrue(candidate["finite_probability_measure_exported"])

    def test_finite_carrier_check(self):
        audit = self.payload["exchangeable_edge_pair_measure_localization_candidate_audit"]
        carrier = audit["finite_carrier_check"]
        self.assertEqual(carrier["decoded_graph_count"], 16828)
        self.assertEqual(carrier["expected_graph_count"], 16828)
        self.assertEqual(carrier["distinct_edge_counts"], [32])
        self.assertEqual(carrier["edge_count_histogram"]["32"], 16828)
        self.assertTrue(carrier["measure_constant_on_full_carrier"])
        self.assertEqual(carrier["constant_edge_count"], 32)

    def test_premise_matrix(self):
        audit = self.payload["exchangeable_edge_pair_measure_localization_candidate_audit"]
        matrix = audit["candidate_premise_matrix"]
        premises = matrix["premises"]
        self.assertTrue(premises["new_typed_object_declared"])
        self.assertTrue(premises["label_gauge_invariant"])
        self.assertFalse(premises["nonconstant_on_full_carrier"])
        self.assertFalse(premises["separates_or_refines_combined_functional"])
        self.assertFalse(premises["field_or_spacetime_support_exported"])
        self.assertTrue(premises["target_independent_units_or_measure_exported"])
        self.assertFalse(premises["coupling_or_variational_chain_rule_exported"])
        self.assertFalse(matrix["accepted_as_strict_localization_or_coupling_object"])

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        facts = matrix["facts"]
        self.assertTrue(facts["new_single_candidate_object_tested"])
        self.assertTrue(facts["candidate_is_label_gauge_invariant"])
        self.assertTrue(facts["finite_probability_measure_exported"])
        self.assertFalse(facts["candidate_nonconstant_on_full_carrier"])
        self.assertFalse(facts["candidate_separates_or_refines_combined_functional"])
        self.assertFalse(facts["field_or_spacetime_support_exported"])
        self.assertFalse(facts["coupling_or_variational_chain_rule_exported"])
        self.assertFalse(facts["selector_bridge_or_role_transfer_imported"])
        self.assertTrue(matrix["accepted_as_finite_label_safe_measure_candidate"])
        self.assertFalse(matrix["accepted_as_strict_localization_or_coupling_object"])
        self.assertTrue(matrix["accepted_as_bounded_candidate_no_go"])

    def test_documents_updated(self):
        self.assertIn("P2842/S1792", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2842/S1792", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2842/S1792", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2842", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
