import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2848_s1798_coupling_coefficient_unit_source_law_audit.py"
JSON_PATH = ROOT / "generated" / "p2848_s1798_coupling_coefficient_unit_source_law_audit.json"
MD_PATH = ROOT / "generated" / "p2848_s1798_coupling_coefficient_unit_source_law_audit.md"


class P2848CouplingCoefficientUnitSourceLawAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status_input_and_carrier(self):
        self.assertEqual(self.payload["status"], "P2848_COUPLING_COEFFICIENT_UNIT_SOURCE_LAW_AUDIT_NO_CLOSURE")
        audit = self.payload["coefficient_unit_source_audit"]
        self.assertEqual(audit["input_statuses_rechecked"]["P2847"], "P2847_TARGET_INDEPENDENT_VOLUME_FORM_UNIT_SOURCE_AUDIT_NO_CLOSURE")
        self.assertTrue(audit["carrier_check"]["coverage_ok"])
        self.assertEqual(audit["carrier_check"]["decoded_graph_count"], 16828)

    def test_mass_stats_and_graph_dependent_normalizers(self):
        stats = self.payload["coefficient_unit_source_audit"]["density_mass_stats"]
        self.assertEqual(stats["uniform_vertex_density"]["distinct_mass_count"], 1)
        self.assertGreater(stats["four_cycle_plus_one_density"]["distinct_mass_count"], 1)
        rows = self.payload["coefficient_unit_source_audit"]["candidate_coefficient_rows"]
        self.assertEqual(rows["uniform_vertex_density"]["inverse_raw_mass"]["distinct_coefficient_count"], 1)
        self.assertGreater(rows["four_cycle_plus_one_density"]["inverse_raw_mass"]["distinct_coefficient_count"], 1)
        self.assertFalse(rows["four_cycle_plus_one_density"]["inverse_raw_mass"]["premises"]["target_independent_across_graphs"])

    def test_no_candidate_exports_unit_source_or_ltotal(self):
        audit = self.payload["coefficient_unit_source_audit"]
        self.assertEqual(audit["accepted_candidate_count"], 0)
        for laws in audit["candidate_coefficient_rows"].values():
            for row in laws.values():
                self.assertFalse(row["premises"]["unit_dimension_source_law"])
                self.assertFalse(row["premises"]["compatible_with_volume_or_pullback"])
                self.assertFalse(row["premises"]["variational_chain_rule"])
                self.assertFalse(row["accepted_as_target_independent_unit_coefficient_source"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["nonproxy_ltotal_term_exported"])
        self.assertFalse(flags["bridge_closure_exported"])
        self.assertFalse(flags["toe_closure_exported"])

    def test_acceptance_and_documents(self):
        self.assertTrue(self.payload["acceptance_matrix"]["accepted_as_coefficient_unit_source_obstruction_audit"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_target_independent_unit_coefficient_source"])
        self.assertIn("P2848/S1798", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2848/S1798", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2848/S1798", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2848", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
