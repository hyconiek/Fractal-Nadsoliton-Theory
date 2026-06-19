import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2913_s1863_gamma_source_provenance_obstruction_audit.py"
JSON_PATH = ROOT / "generated" / "p2913_s1863_gamma_source_provenance_obstruction_audit.json"
MD_PATH = ROOT / "generated" / "p2913_s1863_gamma_source_provenance_obstruction_audit.md"


class P2913GammaSourceProvenanceObstructionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.objects = cls.payload["constructed_theoretical_objects"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_p2912_input(self):
        self.assertEqual(self.payload["status"], "P2913_GAMMA_SOURCE_PROVENANCE_OBSTRUCTION_AUDIT_NO_EXPORT")
        self.assertTrue(self.acceptance["p2912_rechecked_variational_readiness"])

    def test_candidate_family_and_schema(self):
        self.assertEqual(self.acceptance["candidate_source_object_count"], 6)
        self.assertEqual(self.acceptance["accepted_candidate_count"], 0)
        self.assertEqual(len(self.objects["candidate_source_objects"]), 6)
        self.assertEqual(self.objects["missing_theorem_name"], "Strict_Gamma_9_5_Action_Unit_Source_Theorem")
        self.assertGreaterEqual(len(self.objects["acceptance_schema"]), 5)

    def test_generated_scan_has_no_positive_source_hit(self):
        self.assertGreater(self.acceptance["generated_json_files_scanned"], 0)
        self.assertGreaterEqual(self.acceptance["gamma_mention_hit_count"], 1)
        self.assertEqual(self.acceptance["positive_gamma_source_provenance_hit_count"], 0)
        self.assertFalse(self.acceptance["strict_gamma_9_5_source_theorem_exported"])
        self.assertFalse(self.acceptance["accepted_as_nonproxy_ltotal_coupling_source"])

    def test_no_ltotal_or_closure_export(self):
        self.assertFalse(self.acceptance["nonzero_unit_bearing_gamma_value_exported"])
        self.assertFalse(self.acceptance["p2911_p2912_coupling_theorem_exported"])
        self.assertFalse(any(self.flags.values()))

    def test_documents_updated(self):
        self.assertIn("P2913/S1863", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2913/S1863", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2913/S1863", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2913", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
