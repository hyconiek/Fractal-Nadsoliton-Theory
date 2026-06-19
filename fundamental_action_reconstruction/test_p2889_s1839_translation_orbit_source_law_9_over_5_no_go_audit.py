import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2889_s1839_translation_orbit_source_law_9_over_5_no_go_audit.py"
JSON_PATH = ROOT / "generated" / "p2889_s1839_translation_orbit_source_law_9_over_5_no_go_audit.json"
MD_PATH = ROOT / "generated" / "p2889_s1839_translation_orbit_source_law_9_over_5_no_go_audit.md"


class P2889TranslationOrbitSourceLawNoGoAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["translation_orbit_source_law_9_over_5_no_go_audit"]
        cls.orbit = cls.audit["orbit_audit"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2889_TRANSLATION_ORBIT_SOURCE_LAW_9_OVER_5_NO_GO_AUDIT_NO_CLOSURE")
        self.assertEqual(self.audit["input_status_rechecked"], "P2888_NON_C12_ORIGIN_LAW_UNIT_COUPLING_9_OVER_5_NO_GO_AUDIT_NO_CLOSURE")
        self.assertTrue(self.facts["p2888_rechecked"])

    def test_reproduces_p2888_target_count(self):
        self.assertEqual(self.orbit["target_triple_count"], 600)
        self.assertTrue(self.facts["p2888_target_count_reproduced"])

    def test_translation_orbit_obstruction(self):
        self.assertEqual(self.orbit["translation_orbit_count"], 50)
        self.assertEqual(self.orbit["translation_orbit_size_histogram"], {"12": 50})
        self.assertEqual(self.orbit["translation_stabilizer_size_histogram"], {"1": 50})
        self.assertTrue(self.orbit["all_target_triples_have_free_translation_orbit"])
        self.assertTrue(self.facts["translation_neutral_rule_can_select_only_orbit_not_representative"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_translation_neutral_strict_source_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_unique_embedded_support_origin_density_triple"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_nonimported_9_over_5_variational_chain_rule"])

    def test_documents_updated(self):
        self.assertIn("P2889/S1839", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2889/S1839", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2889/S1839", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2889", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
