import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2847_s1797_target_independent_volume_form_unit_source_audit.py"
JSON_PATH = ROOT / "generated" / "p2847_s1797_target_independent_volume_form_unit_source_audit.json"
MD_PATH = ROOT / "generated" / "p2847_s1797_target_independent_volume_form_unit_source_audit.md"


class P2847TargetIndependentVolumeFormUnitSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2847_TARGET_INDEPENDENT_VOLUME_FORM_UNIT_SOURCE_AUDIT_NO_CLOSURE")
        audit = self.payload["volume_form_unit_source_audit"]
        self.assertEqual(audit["input_statuses_rechecked"]["P2846"], "P2846_LABEL_SAFE_VERTEX_LOCALIZATION_PULLBACK_CANDIDATE_NO_GO_NO_CLOSURE")
        self.assertTrue(audit["carrier_check"]["coverage_ok"])
        self.assertEqual(audit["carrier_check"]["decoded_graph_count"], 16828)

    def test_candidate_volume_stats(self):
        rows = self.payload["volume_form_unit_source_audit"]["candidate_rows"]
        self.assertEqual(len(rows), 4)
        self.assertTrue(rows["uniform_vertex_density"]["premises"]["target_independent_total_volume"])
        self.assertFalse(rows["uniform_vertex_density"]["premises"]["unit_dimension_source"])
        self.assertGreater(rows["four_cycle_plus_one_density"]["finite_stats"]["distinct_raw_total_mass_count"], 1)
        self.assertFalse(rows["four_cycle_plus_one_density"]["premises"]["target_independent_total_volume"])
        self.assertGreater(rows["triangle_square_plus_one_density"]["finite_stats"]["distinct_normalized_probability_profile_count"], 1)

    def test_no_candidate_exports_units_or_ltotal(self):
        audit = self.payload["volume_form_unit_source_audit"]
        self.assertEqual(audit["accepted_candidate_count"], 0)
        for row in audit["candidate_rows"].values():
            self.assertFalse(row["premises"]["unit_dimension_source"])
            self.assertFalse(row["premises"]["canonical_vertex_to_field_support"])
            self.assertFalse(row["premises"]["coupling_coefficient_rule"])
            self.assertFalse(row["accepted_as_target_independent_unit_volume_source"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["nonproxy_ltotal_term_exported"])
        self.assertFalse(flags["eom_closure_exported"])
        self.assertFalse(flags["hamiltonian_closure_exported"])

    def test_acceptance_and_documents(self):
        self.assertTrue(self.payload["acceptance_matrix"]["accepted_as_volume_form_unit_source_obstruction_audit"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_target_independent_unit_volume_source"])
        self.assertIn("P2847/S1797", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2847/S1797", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2847/S1797", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2847", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
