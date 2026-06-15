import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2779_s1729_16node_circulant_full_spectrum_quotient_audit.py"
OUT = ROOT / "generated" / "p2779_s1729_16node_circulant_full_spectrum_quotient_audit.json"
MD = ROOT / "generated" / "p2779_s1729_16node_circulant_full_spectrum_quotient_audit.md"


class P2779FullSpectrum16NodeQuotientAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2779_16NODE_CIRCULANT_FULL_SPECTRUM_QUOTIENT_AUDIT_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2778"], "P2778_MAX_SYMMETRY_16NODE_GEOMETRY_SOURCE_AUDIT_NO_CLOSURE")
        self.assertIn("full Laplacian spectrum", self.payload["audited_question"])

    def test_quotient_full_spectrum_has_no_nonisomorphic_collision(self):
        witness = self.payload["full_spectrum_16node_quotient_witness"]
        self.assertEqual(witness["labeled_candidate_count"], 19)
        self.assertEqual(witness["isomorphism_class_count"], 6)
        self.assertEqual(witness["distinct_laplacian_spectrum_count_after_quotient"], 6)
        self.assertEqual(witness["nonisomorphic_cospectral_collision_count"], 0)
        self.assertTrue(witness["full_spectrum_injective_on_declared_quotient"])
        self.assertTrue(witness["torus_4x4_representative_present"])
        representatives = {row["representative"] for row in witness["quotient_rows"]}
        self.assertIn("torus_4x4", representatives)
        self.assertIn("circulant_pm1_pm7", representatives)

    def test_acceptance_is_positive_but_blocks_closure(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["facts"]["declared_16node_class_audited"])
        self.assertTrue(acceptance["facts"]["isomorphism_quotient_performed"])
        self.assertTrue(acceptance["accepted_as_declared_class_full_spectrum_uniqueness_witness"])
        self.assertFalse(acceptance["accepted_as_strict_spectral_source_law"])
        self.assertFalse(acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("strict_nadsoliton_spectral_source_law_exported", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("broader 16-node 4-regular non-circulant", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2779", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2779/S1729", MD.read_text(encoding="utf-8"))
        self.assertIn("P2779/S1729", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2779/S1729", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2779", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
