import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2850_s1800_eml_single_operator_kernel_bridge_impact_audit.py"
JSON_PATH = ROOT / "generated" / "p2850_s1800_eml_single_operator_kernel_bridge_impact_audit.json"
MD_PATH = ROOT / "generated" / "p2850_s1800_eml_single_operator_kernel_bridge_impact_audit.md"


class P2850EmlSingleOperatorKernelBridgeImpactAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status_input_and_sources(self):
        self.assertEqual(self.payload["status"], "P2850_EML_SINGLE_OPERATOR_KERNEL_BRIDGE_IMPACT_AUDIT_NO_CLOSURE")
        audit = self.payload["eml_impact_audit"]
        self.assertEqual(audit["input_statuses_rechecked"]["P2849"], "P2849_DAMPING_COMPRESSION_KERNEL_BRIDGE_ATOM_AUDIT_NO_CLOSURE")
        sources = self.payload["external_sources_checked"]
        self.assertEqual(sources["arxiv_id_version"], "2603.21852v2")
        self.assertIn("arxiv.org", sources["arxiv_html"])

    def test_eml_syntax_positive_but_bridge_negative(self):
        audit = self.payload["eml_impact_audit"]
        self.assertLess(audit["max_exp_identity_abs_error"], 1e-15)
        classification = audit["elementary_syntax_classification"]
        self.assertTrue(classification["legacy_formula_elementary"])
        self.assertTrue(classification["strict_formula_elementary"])
        self.assertEqual(classification["eml_representation_relevance"], "syntax_only_unless_typed_source_law_added")
        self.assertGreater(audit["max_legacy_strict_kernel_abs_gap_on_audited_distances"], 0.0)

    def test_missing_source_premises_and_no_closure(self):
        matrix = self.payload["eml_impact_audit"]["premise_matrix"]
        self.assertFalse(matrix["changes_p2849_bridge_boundary"])
        self.assertIn("eml_exports_beta_eta_source", matrix["missing_premises"])
        self.assertIn("eml_exports_completion_map_semantics", matrix["missing_premises"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["kernel_bridge_exported"])
        self.assertFalse(flags["beta_eta_source_exported"])
        self.assertFalse(flags["ltotal_exported"])
        self.assertFalse(flags["toe_closure_exported"])

    def test_acceptance_and_documents(self):
        self.assertTrue(self.payload["acceptance_matrix"]["accepted_as_eml_external_information_impact_audit"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_new_kernel_bridge_source"])
        self.assertIn("P2850/S1800", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2850/S1800", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2850/S1800", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2850", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
