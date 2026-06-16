import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
JSON_PATH = ROOT / "generated" / "p2800_s1750_girth4_shortcode_import_absence_manifest.json"
MD_PATH = ROOT / "generated" / "p2800_s1750_girth4_shortcode_import_absence_manifest.md"


class P2800Girth4ShortcodeImportAbsenceManifestTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON_PATH.exists():
            import subprocess, sys
            subprocess.run([sys.executable, str(ROOT / "p2800_s1750_girth4_shortcode_import_absence_manifest.py")], check=True)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.manifest = cls.payload["girth4_shortcode_import_absence_manifest"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2800_GIRTH4_SHORTCODE_IMPORT_ABSENCE_MANIFEST_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2798"], "P2798_EXTERNAL_GIRTH4_SUBTARGET_LOCAL_GIRTH_FILTER_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2799"], "P2799_LOCAL_GIRTH4_SPECTRAL_ORBIT_TABLE_CERTIFICATE_NO_CLOSURE")

    def test_expected_target_and_absence_scan(self):
        self.assertEqual(self.manifest["expected_girth4_shortcode_class_count"], 16828)
        self.assertEqual(self.manifest["p2798_external_target_count"], 16828)
        self.assertEqual(self.manifest["p2799_local_girth4_table_row_count"], 6)
        self.assertEqual(self.manifest["exact_16828_line_candidate_file_count"], 0)
        self.assertFalse(self.manifest["required_shortcode_artifact_present"])
        self.assertEqual(self.manifest["subtarget_gap_if_no_import"], 16822)

    def test_candidate_scan_has_generated_root(self):
        roots = {row["root"] for row in self.manifest["candidate_scan_rows"]}
        self.assertIn("fundamental_action_reconstruction/generated", roots)
        self.assertIn(".json", self.manifest["candidate_suffixes"])

    def test_acceptance_blocks_closure(self):
        self.assertTrue(self.acceptance["accepted_as_girth4_shortcode_import_absence_manifest"])
        self.assertFalse(self.acceptance["accepted_as_girth4_shortcode_graph_list_import"])
        self.assertFalse(self.acceptance["accepted_as_full_16node_canonical_generator_certificate"])
        self.assertFalse(self.acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(self.acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("required_shortcode_artifact_present", self.acceptance["missing_criteria"])
        self.assertIn("P2697-P2800", self.payload["decision"]["next_honest_step"])

    def test_markdown_and_guardrails_written(self):
        self.assertIn("P2800/S1750", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2800/S1750", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2800/S1750", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2800", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
