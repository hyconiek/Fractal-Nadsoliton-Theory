import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2891_s1841_strict_phase_origin_source_artifact_inventory_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2891_s1841_strict_phase_origin_source_artifact_inventory_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2891_s1841_strict_phase_origin_source_artifact_inventory_no_go_audit.md"


class P2891StrictPhaseOriginSourceArtifactInventoryNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["strict_phase_origin_source_artifact_inventory_no_go_audit"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2891_STRICT_PHASE_ORIGIN_SOURCE_ARTIFACT_INVENTORY_NO_GO_AUDIT_NO_CLOSURE")
        self.assertEqual(self.audit["input_status_rechecked"], "P2890_FOURIER_PHASE_SOURCE_LAW_9_OVER_5_NO_GO_AUDIT_NO_CLOSURE")
        self.assertTrue(self.facts["p2890_rechecked"])

    def test_inventory_counts(self):
        self.assertGreaterEqual(self.audit["generated_json_file_count_excluding_self"], 4938)
        self.assertGreater(self.audit["relevant_record_count"], 0)
        self.assertEqual(self.audit["coupled_positive_export_record_count"], 0)
        self.assertTrue(self.facts["generated_json_inventory_nonempty"])
        self.assertTrue(self.facts["phase_or_coupling_terms_found"])
        self.assertTrue(self.facts["no_coupled_positive_phase_origin_9_over_5_export_found"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_strict_phase_origin_source_artifact"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_nonconventional_phase_or_sign_source"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_coupling_to_9_over_5_variational_density"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_nonimported_9_over_5_variational_chain_rule"])

    def test_documents_updated(self):
        self.assertIn("P2891/S1841", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2891/S1841", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2891/S1841", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2891", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
