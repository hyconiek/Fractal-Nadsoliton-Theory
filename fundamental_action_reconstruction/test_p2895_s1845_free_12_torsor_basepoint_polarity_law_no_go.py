import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2895_s1845_free_12_torsor_basepoint_polarity_law_no_go.py"
JSON_PATH = ROOT / "generated" / "p2895_s1845_free_12_torsor_basepoint_polarity_law_no_go.json"
MD_PATH = ROOT / "generated" / "p2895_s1845_free_12_torsor_basepoint_polarity_law_no_go.md"


class P2895Free12TorsorBasepointPolarityLawNoGoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["free_12_torsor_basepoint_polarity_law_no_go"]["torsor_audit"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2895_FREE_12_TORSOR_BASEPOINT_POLARITY_LAW_NO_GO_NO_CLOSURE")
        self.assertTrue(self.facts["p2894_rechecked"])
        self.assertEqual(self.payload["free_12_torsor_basepoint_polarity_law_no_go"]["input_status_rechecked"], "P2894_TRANSLATION_BREAKING_SOURCE_ALPHABET_CAPACITY_NO_GO_NO_CLOSURE")

    def test_torsor_counts(self):
        self.assertEqual(self.audit["source_torsor_size"], 12)
        self.assertEqual(self.audit["target_torsor_size"], 12)
        self.assertEqual(self.audit["equivariant_map_count"], 12)
        self.assertEqual(self.audit["equivariant_map_offsets"], list(range(12)))
        self.assertEqual(self.audit["invariant_basepoint_count"], 0)
        self.assertEqual(self.audit["self_opposite_offsets"], [0, 6])
        self.assertEqual(self.audit["nonself_opposite_pair_count"], 5)

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_strict_free_12_torsor_source_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_nonimported_basepoint_or_polarity_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_unique_coupling_to_9_over_5_carrier"])

    def test_documents_updated(self):
        self.assertIn("P2895/S1845", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2895/S1845", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2895/S1845", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2895", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
