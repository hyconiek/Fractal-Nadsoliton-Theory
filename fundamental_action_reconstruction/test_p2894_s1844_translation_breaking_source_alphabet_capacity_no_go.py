import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2894_s1844_translation_breaking_source_alphabet_capacity_no_go.py"
JSON_PATH = ROOT / "generated" / "p2894_s1844_translation_breaking_source_alphabet_capacity_no_go.json"
MD_PATH = ROOT / "generated" / "p2894_s1844_translation_breaking_source_alphabet_capacity_no_go.md"


class P2894TranslationBreakingSourceAlphabetCapacityNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["translation_breaking_source_alphabet_capacity_no_go"]["capacity_audit"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2894_TRANSLATION_BREAKING_SOURCE_ALPHABET_CAPACITY_NO_GO_NO_CLOSURE")
        self.assertTrue(self.facts["p2893_rechecked"])
        self.assertEqual(self.payload["translation_breaking_source_alphabet_capacity_no_go"]["input_status_rechecked"], "P2893_FREE_ORBIT_SECTION_OBSTRUCTION_THEOREM_NO_CLOSURE")

    def test_low_capacity_exhaustion(self):
        self.assertEqual(self.audit["low_capacity_source_alphabet_sizes_tested"], list(range(1, 12)))
        self.assertGreater(self.audit["low_capacity_z12_set_type_count"], 0)
        self.assertEqual(self.audit["low_capacity_types_with_equivariant_map_to_free_12_orbit"], 0)
        self.assertEqual(self.audit["minimum_source_orbit_size_required"], 12)
        self.assertTrue(self.audit["free_12_torsor_necessary_for_equivariant_targeting"])

    def test_named_source_rows(self):
        rows = {row["name"]: row for row in self.audit["named_source_capacity_rows"]}
        for name in ["scalar/trivial", "binary sign", "ternary phase", "quarter phase", "half-cycle phase"]:
            self.assertEqual(rows[name]["equivariant_map_count_to_free_12_orbit"], 0)
        self.assertEqual(rows["full Z12 phase torsor"]["equivariant_map_count_to_free_12_orbit"], 12)

    def test_false_export_flags_and_documents(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_translation_breaking_strict_source"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_basepoint_or_polarity_law"])
        self.assertIn("P2894/S1844", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2894/S1844", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2894/S1844", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2894", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
