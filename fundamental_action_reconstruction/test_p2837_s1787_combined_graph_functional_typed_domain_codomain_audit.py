import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2837_s1787_combined_graph_functional_typed_domain_codomain_audit.py"
JSON_PATH = ROOT / "generated" / "p2837_s1787_combined_graph_functional_typed_domain_codomain_audit.json"
MD_PATH = ROOT / "generated" / "p2837_s1787_combined_graph_functional_typed_domain_codomain_audit.md"


class P2837CombinedGraphFunctionalTypedDomainCodomainAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status(self):
        self.assertEqual(
            self.payload["status"],
            "P2837_TYPED_DOMAIN_CODOMAIN_OBLIGATION_NO_GO_NO_COUPLING_NO_CLOSURE",
        )

    def test_separator_rechecked(self):
        audit = self.payload["combined_graph_functional_typed_domain_codomain_audit"]
        sep = audit["combined_separator_rechecked"]
        self.assertEqual(sep["combined_class_count"], 16828)
        self.assertEqual(sep["combined_collision_class_count"], 0)
        self.assertEqual(sep["p2834_patch_graph_count"], 138)

    def test_candidate_typed_maps(self):
        audit = self.payload["combined_graph_functional_typed_domain_codomain_audit"]
        candidates = {row["candidate"]: row for row in audit["candidate_typed_maps"]}
        self.assertEqual(set(candidates), {
            "graph_to_kernel_scalar_control",
            "graph_to_lagrangian_density_term",
            "graph_to_source_density",
            "graph_to_dimensionless_coefficient",
        })
        self.assertTrue(all(row["premises"]["source_domain_object"] for row in candidates.values()))
        self.assertTrue(all(row["premises"]["target_codomain_object"] for row in candidates.values()))
        self.assertFalse(any(row["accepted_as_typed_K_or_Ltotal_map"] for row in candidates.values()))
        self.assertTrue(candidates["graph_to_dimensionless_coefficient"]["premises"]["evaluation_or_pullback_map"])
        self.assertFalse(candidates["graph_to_dimensionless_coefficient"]["premises"]["target_independent_units_or_scale"])

    def test_obstruction_summary(self):
        audit = self.payload["combined_graph_functional_typed_domain_codomain_audit"]
        summary = audit["typed_obstruction_summary"]
        self.assertEqual(summary["candidate_count"], 4)
        self.assertEqual(summary["accepted_typed_map_count"], 0)
        self.assertEqual(summary["accepted_typed_maps"], [])
        self.assertIn("target_independent_units_or_scale", summary["common_hard_blockers"])
        self.assertIn("locality_or_covariance_rule", summary["common_hard_blockers"])
        self.assertIn("variational_variable_identification", summary["common_hard_blockers"])
        self.assertIn("coupling_coefficient_rule", summary["common_hard_blockers"])
        self.assertEqual(summary["blocker_histogram"]["evaluation_or_pullback_map"], 3)

    def test_acceptance_boundaries(self):
        matrix = self.payload["acceptance_matrix"]
        facts = matrix["facts"]
        self.assertTrue(facts["p2835_combined_separator_rechecked"])
        self.assertTrue(facts["at_least_one_candidate_typed_map_tested"])
        self.assertFalse(facts["accepted_typed_K_or_Ltotal_map_exported"])
        self.assertFalse(facts["target_independent_units_available_from_p2836"])
        self.assertFalse(facts["selector_bridge_or_role_transfer_imported"])
        self.assertTrue(matrix["accepted_as_typed_domain_codomain_audit"])
        self.assertFalse(matrix["accepted_as_typed_K_or_Ltotal_map"])
        self.assertTrue(matrix["accepted_as_typed_domain_codomain_no_go"])

    def test_documents_updated(self):
        self.assertIn("P2837/S1787", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2837/S1787", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2837/S1787", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2837", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
