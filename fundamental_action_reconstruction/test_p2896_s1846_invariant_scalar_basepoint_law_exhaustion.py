import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2896_s1846_invariant_scalar_basepoint_law_exhaustion.py"
JSON_PATH = ROOT / "generated" / "p2896_s1846_invariant_scalar_basepoint_law_exhaustion.json"
MD_PATH = ROOT / "generated" / "p2896_s1846_invariant_scalar_basepoint_law_exhaustion.md"


class P2896InvariantScalarBasepointLawExhaustionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["invariant_scalar_basepoint_law_exhaustion"]["scalar_exhaustion"]
        cls.facts = cls.payload["acceptance_matrix"]["facts"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2896_INVARIANT_SCALAR_BASEPOINT_LAW_EXHAUSTION_NO_CLOSURE")
        self.assertTrue(self.facts["p2895_rechecked"])
        self.assertEqual(
            self.payload["invariant_scalar_basepoint_law_exhaustion"]["input_status_rechecked"],
            "P2895_FREE_12_TORSOR_BASEPOINT_POLARITY_LAW_NO_GO_NO_CLOSURE",
        )

    def test_scalar_exhaustion_counts(self):
        self.assertEqual(self.audit["torsor_size"], 12)
        self.assertEqual(self.audit["tested_trivial_alphabet_sizes"], [1, 2, 3, 4, 6, 12])
        self.assertEqual(self.audit["total_invariant_scalar_law_count"], 28)
        self.assertEqual(self.audit["total_unique_marked_level_law_count"], 0)
        self.assertEqual(self.audit["total_unique_argmin_law_count"], 0)
        self.assertEqual(self.audit["total_unique_argmax_law_count"], 0)
        self.assertTrue(all(row["all_invariant_laws_constant"] for row in self.audit["alphabet_rows"]))

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.payload["acceptance_matrix"]["accepted_as_nonimported_basepoint_law"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_translation_breaking_scalar_source"])
        self.assertFalse(self.payload["acceptance_matrix"]["exports_unique_coupling_to_9_over_5_carrier"])

    def test_documents_updated(self):
        self.assertIn("P2896/S1846", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2896/S1846", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2896/S1846", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2896", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
