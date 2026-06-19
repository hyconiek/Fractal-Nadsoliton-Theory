import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2908_s1858_joint_source_provenance_alternative_audit.py"
JSON_PATH = ROOT / "generated" / "p2908_s1858_joint_source_provenance_alternative_audit.json"
MD_PATH = ROOT / "generated" / "p2908_s1858_joint_source_provenance_alternative_audit.md"


class P2908JointSourceProvenanceAlternativeAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.objects = cls.payload["constructed_theoretical_objects"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2908_JOINT_SOURCE_PROVENANCE_ALTERNATIVE_AUDIT_NO_STRICT_PROVENANCE")
        self.assertTrue(self.acceptance["p2907_rechecked_candidate_not_strict"])
        self.assertTrue(self.acceptance["distinguished_j_0_plus_present"])

    def test_alternative_family_counts(self):
        self.assertEqual(self.acceptance["alternative_count"], 24)
        self.assertEqual(self.acceptance["translated_sign_flipped_alternative_count"], 23)
        self.assertEqual(len(self.objects["joint_source_alternative_family"]), 24)
        j0 = self.objects["distinguished_candidate_rechecked"]
        self.assertEqual(j0["origin"], 0)
        self.assertEqual(j0["sign"], 1)
        self.assertEqual(j0["defect_edge"], [0, 5])

    def test_no_positive_joint_provenance_hit(self):
        scan = self.objects["provenance_scan"]
        self.assertGreater(scan["generated_json_file_count"], 0)
        self.assertGreaterEqual(scan["mention_hit_count"], 1)
        self.assertEqual(scan["positive_provenance_hit_count"], 0)
        self.assertEqual(self.acceptance["positive_provenance_hit_count"], 0)
        self.assertFalse(self.acceptance["strict_provenance_for_j_0_plus_exported"])
        self.assertFalse(self.acceptance["accepted_as_strict_joint_source"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.flags["j_0_plus_provenance_exported"])
        self.assertFalse(self.flags["nonproxy_ltotal_exported"])
        self.assertFalse(self.flags["toe_closure_exported"])

    def test_documents_updated(self):
        self.assertIn("P2908/S1858", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2908/S1858", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2908/S1858", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2908", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
