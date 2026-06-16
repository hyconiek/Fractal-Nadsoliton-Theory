import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GEN = ROOT / "generated"
JSON = GEN / "p2785_s1735_exact_characteristic_polynomial_certificate.json"
MD = GEN / "p2785_s1735_exact_characteristic_polynomial_certificate.md"


class P2785ExactCharacteristicPolynomialCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not JSON.exists():
            import p2785_s1735_exact_characteristic_polynomial_certificate as p2785
            p2785.main()
        cls.payload = json.loads(JSON.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2785_EXACT_CHARACTERISTIC_POLYNOMIAL_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2784"], "P2784_SEVEN_CLASS_MULTISPECTRAL_ROBUSTNESS_AUDIT_NO_CLOSURE")
        self.assertIn("exact integer characteristic-polynomial", self.payload["audited_question"])

    def test_exact_polynomial_counts(self):
        witness = self.payload["exact_polynomial_witness"]
        self.assertEqual(witness["representative_count"], 7)
        self.assertEqual(witness["pair_count"], 21)
        self.assertEqual(witness["exact_charpoly_collision_counts"], {
            "adjacency_charpoly_coefficients": 0,
            "laplacian_charpoly_coefficients": 0,
            "signless_laplacian_charpoly_coefficients": 0,
        })
        self.assertTrue(witness["all_pairs_separated_by_all_three_exact_charpolys"])
        self.assertEqual(len(witness["exact_rows"]), 7)
        self.assertEqual(len(witness["pair_rows"]), 21)
        for row in witness["exact_rows"]:
            self.assertEqual(len(row["adjacency_charpoly_coefficients"]), 17)
            self.assertEqual(len(row["laplacian_charpoly_coefficients"]), 17)
            self.assertEqual(len(row["signless_laplacian_charpoly_coefficients"]), 17)
            self.assertGreater(row["kirchhoff_spanning_tree_count_exact"], 0)
        for row in witness["pair_rows"]:
            self.assertTrue(row["all_three_charpolys_distinct"])

    def test_acceptance_blocks_closure(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["accepted_as_exact_local_algebra_certificate"])
        self.assertFalse(acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        for flag, value in self.payload["decision"]["negative_export_flags"].items():
            self.assertFalse(value, flag)
        self.assertIn("P2697-P2785", self.payload["decision"]["next_honest_step"])

    def test_documentation_updated(self):
        self.assertIn("P2785/S1735", MD.read_text(encoding="utf-8"))
        self.assertIn("P2785/S1735", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2785/S1735", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2785", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
