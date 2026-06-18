import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2838_s1788_combined_graph_functional_variational_derivative_obstruction_audit.py"
JSON_PATH = ROOT / "generated" / "p2838_s1788_combined_graph_functional_variational_derivative_obstruction_audit.json"
MD_PATH = ROOT / "generated" / "p2838_s1788_combined_graph_functional_variational_derivative_obstruction_audit.md"


class P2838CombinedGraphFunctionalVariationalDerivativeObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status(self):
        self.assertEqual(
            self.payload["status"],
            "P2838_VARIATIONAL_DERIVATIVE_OBSTRUCTION_NO_GO_NO_COUPLING_NO_CLOSURE",
        )

    def test_finite_difference_inventory(self):
        audit = self.payload["combined_graph_functional_variational_derivative_obstruction_audit"]
        inventory = audit["finite_difference_inventory"]
        self.assertEqual(inventory["adjacency_bit_variable_count"], 120)
        self.assertEqual(inventory["first_difference_slots_per_graph"], 120)
        self.assertEqual(inventory["second_difference_slots_per_graph"], 7140)
        self.assertTrue(inventory["p2833_first_variation_witness_available"])
        self.assertTrue(inventory["p2834_second_variation_witness_available_on_residuals"])

    def test_candidate_variational_routes(self):
        audit = self.payload["combined_graph_functional_variational_derivative_obstruction_audit"]
        routes = {row["candidate"]: row for row in audit["candidate_variational_derivative_routes"]}
        self.assertEqual(set(routes), {
            "adjacency_bit_discrete_gradient",
            "field_variation_delta_F_delta_phi",
            "metric_variation_delta_F_delta_g",
            "kernel_parameter_variation_delta_F_delta_K",
        })
        self.assertTrue(all(row["premises"]["finite_difference_operator_on_graph_bits"] for row in routes.values()))
        self.assertFalse(any(row["accepted_as_formal_variational_derivative"] for row in routes.values()))
        self.assertFalse(routes["adjacency_bit_discrete_gradient"]["premises"]["continuous_or_formal_field_variable"])
        self.assertTrue(routes["field_variation_delta_F_delta_phi"]["premises"]["continuous_or_formal_field_variable"])

    def test_obstruction_summary(self):
        audit = self.payload["combined_graph_functional_variational_derivative_obstruction_audit"]
        summary = audit["variational_obstruction_summary"]
        self.assertEqual(summary["candidate_route_count"], 4)
        self.assertEqual(summary["accepted_formal_variational_derivative_count"], 0)
        self.assertEqual(summary["accepted_routes"], [])
        self.assertIn("localization_or_pullback_Aij_to_fields", summary["common_hard_blockers"])
        self.assertIn("action_integral_or_density_embedding", summary["common_hard_blockers"])
        self.assertIn("target_independent_units_or_coupling_coefficient", summary["common_hard_blockers"])
        self.assertIn("chain_rule_from_graph_bits_to_K_or_Ltotal", summary["common_hard_blockers"])
        self.assertEqual(summary["blocker_histogram"]["continuous_or_formal_field_variable"], 1)

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        facts = matrix["facts"]
        self.assertTrue(facts["p2835_combined_separator_rechecked"])
        self.assertTrue(facts["finite_graph_difference_operators_available"])
        self.assertFalse(facts["accepted_formal_variational_derivative_exported"])
        self.assertFalse(facts["typed_domain_codomain_available_from_p2837"])
        self.assertFalse(facts["selector_bridge_or_role_transfer_imported"])
        self.assertTrue(matrix["accepted_as_finite_difference_inventory"])
        self.assertFalse(matrix["accepted_as_formal_variational_derivative"])
        self.assertTrue(matrix["accepted_as_variational_derivative_no_go"])

    def test_documents_updated(self):
        self.assertIn("P2838/S1788", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2838/S1788", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2838/S1788", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2838", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
