import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2783_s1733_seven_class_quotient_integrity_certificate.py"
OUT = ROOT / "generated" / "p2783_s1733_seven_class_quotient_integrity_certificate.json"
MD = ROOT / "generated" / "p2783_s1733_seven_class_quotient_integrity_certificate.md"


class P2783SevenClassQuotientIntegrityCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2783_SEVEN_CLASS_QUOTIENT_INTEGRITY_CERTIFICATE_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2781"], "P2781_ENUMERATED_TWO_LAYER_C8_SPECTRUM_COLLISION_AUDIT_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2782"], "P2782_BIPARTITE_REGULAR_ENUMERATOR_SCALE_OBSTRUCTION_NO_CLOSURE")
        self.assertIn("seven P2781 quotient classes", self.payload["audited_question"])

    def test_all_seven_classes_are_pairwise_clean(self):
        witness = self.payload["seven_class_integrity_witness"]
        self.assertEqual(witness["representative_count"], 7)
        self.assertEqual(witness["pair_count"], 21)
        self.assertEqual(witness["direct_isomorphism_collision_count"], 0)
        self.assertEqual(witness["full_spectrum_collision_count"], 0)
        self.assertTrue(witness["all_representatives_pairwise_nonisomorphic"])
        self.assertTrue(witness["all_representatives_pairwise_spectrally_distinct"])
        self.assertTrue(all(row["separation_certificate"] == "direct_nonisomorphism_and_distinct_spectrum" for row in witness["pair_rows"]))

    def test_acceptance_is_integrity_certificate_not_closure(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["facts"]["all_21_pairs_checked_by_direct_isomorphism"])
        self.assertTrue(acceptance["facts"]["all_representatives_pairwise_nonisomorphic"])
        self.assertTrue(acceptance["facts"]["all_representatives_pairwise_spectrally_distinct"])
        self.assertTrue(acceptance["accepted_as_p2781_quotient_integrity_certificate"])
        self.assertFalse(acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("strict_nadsoliton_spectral_source_law_exported", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("local quotient-integrity certificate", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2783", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2783/S1733", MD.read_text(encoding="utf-8"))
        self.assertIn("P2783/S1733", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2783/S1733", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2783", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
