import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
JSON = GEN / "p2784_s1734_seven_class_multispectral_robustness_audit.json"
MD = GEN / "p2784_s1734_seven_class_multispectral_robustness_audit.md"


class P2784SevenClassMultispectralRobustnessAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON.exists():
            import p2784_s1734_seven_class_multispectral_robustness_audit as p2784
            p2784.main()
        cls.payload = json.loads(JSON.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2784_SEVEN_CLASS_MULTISPECTRAL_ROBUSTNESS_AUDIT_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2783"], "P2783_SEVEN_CLASS_QUOTIENT_INTEGRITY_CERTIFICATE_NO_CLOSURE")
        self.assertIn("adjacency, Laplacian, and signless-Laplacian", self.payload["audited_question"])

    def test_multispectral_counts(self):
        witness = self.payload["multispectral_witness"]
        self.assertEqual(witness["representative_count"], 7)
        self.assertEqual(witness["pair_count"], 21)
        self.assertEqual(witness["spectral_collision_counts"], {
            "adjacency_spectrum": 0,
            "laplacian_spectrum": 0,
            "signless_laplacian_spectrum": 0,
        })
        self.assertTrue(witness["all_pairs_separated_by_all_three_spectra"])
        self.assertEqual(len(witness["invariant_rows"]), 7)
        self.assertEqual(len(witness["pair_rows"]), 21)
        for row in witness["pair_rows"]:
            self.assertTrue(row["all_three_spectra_distinct"])

    def test_acceptance_blocks_closure(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["accepted_as_multispectral_robustness_certificate"])
        self.assertFalse(acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        for flag, value in self.payload["decision"]["negative_export_flags"].items():
            self.assertFalse(value, flag)
        self.assertIn("P2697-P2784", self.payload["decision"]["next_honest_step"])

    def test_documentation_updated(self):
        self.assertIn("P2784/S1734", MD.read_text(encoding="utf-8"))
        self.assertIn("P2784/S1734", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2784/S1734", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2784", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
