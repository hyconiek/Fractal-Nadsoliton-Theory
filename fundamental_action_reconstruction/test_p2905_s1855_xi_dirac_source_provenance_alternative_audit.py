import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2905_s1855_xi_dirac_source_provenance_alternative_audit.py"
JSON_PATH = ROOT / "generated" / "p2905_s1855_xi_dirac_source_provenance_alternative_audit.json"
MD_PATH = ROOT / "generated" / "p2905_s1855_xi_dirac_source_provenance_alternative_audit.md"


class P2905XiDiracSourceProvenanceAlternativeAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.objects = cls.payload["constructed_theoretical_objects"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2905_XI_DIRAC_SOURCE_PROVENANCE_ALTERNATIVE_AUDIT_NO_STRICT_PROVENANCE")
        self.assertTrue(self.acceptance["p2904_rechecked"])
        self.assertTrue(self.acceptance["candidate_xi_0_plus_present"])

    def test_alternative_family_counts(self):
        self.assertEqual(self.acceptance["alternative_count"], 24)
        self.assertEqual(self.acceptance["translated_sign_flipped_alternative_count"], 23)
        self.assertEqual(len(self.objects["alternative_family"]), 24)
        xi = self.objects["distinguished_candidate_rechecked"]
        self.assertEqual(xi["basepoint"], 0)
        self.assertEqual(xi["sign"], 1)
        self.assertEqual(xi["defect_edge"], [0, 5])

    def test_no_positive_provenance_hit(self):
        scan = self.objects["provenance_scan"]
        self.assertGreater(scan["generated_json_file_count"], 0)
        self.assertGreaterEqual(scan["mention_hit_count"], 1)
        self.assertEqual(scan["positive_provenance_hit_count"], 0)
        self.assertEqual(self.acceptance["positive_provenance_hit_count"], 0)
        self.assertFalse(self.acceptance["strict_provenance_for_xi_0_plus_exported"])
        self.assertFalse(self.acceptance["accepted_as_strict_source"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.flags["xi_0_plus_provenance_exported"])
        self.assertFalse(self.flags["nonproxy_ltotal_exported"])
        self.assertFalse(self.flags["toe_closure_exported"])

    def test_documents_updated(self):
        self.assertIn("P2905/S1855", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2905/S1855", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2905/S1855", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2905", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
