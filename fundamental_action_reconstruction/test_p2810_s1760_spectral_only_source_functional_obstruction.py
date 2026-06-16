import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
JSON_PATH = ROOT / "generated" / "p2810_s1760_spectral_only_source_functional_obstruction.json"
MD_PATH = ROOT / "generated" / "p2810_s1760_spectral_only_source_functional_obstruction.md"


class P2810SpectralOnlySourceFunctionalObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(ROOT / "p2810_s1760_spectral_only_source_functional_obstruction.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.matrix = cls.payload["spectral_only_obstruction_matrix"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2810_SPECTRAL_PAIR_ONLY_SOURCE_FUNCTIONAL_OBSTRUCTION_NO_STRICT_SOURCE_LAW_NO_CLOSURE")
        self.assertIn("P2804", self.payload["input_statuses"])
        self.assertIn("P2805", self.payload["input_statuses"])
        self.assertIn("P2808", self.payload["input_statuses"])
        self.assertIn("P2809", self.payload["input_statuses"])

    def test_exact_obstruction_counts(self):
        self.assertEqual(self.matrix["canonical_certificate_class_count_from_p2808"], 16828)
        self.assertEqual(self.matrix["spectral_pair_class_count_from_p2804"], 16211)
        self.assertEqual(self.matrix["quotient_defect_canonical_minus_spectral_pair"], 617)
        self.assertEqual(self.matrix["spectral_pair_collision_class_count"], 578)
        self.assertEqual(self.matrix["spectral_pair_collision_graph_count"], 1195)
        self.assertEqual(self.matrix["p2805_positive_isomorphism_matches_inside_collisions"], 0)
        self.assertEqual(self.matrix["p2805_negative_isomorphism_rejections_inside_collisions"], 660)

    def test_acceptance_boundaries(self):
        self.assertTrue(self.acceptance["accepted_as_spectral_pair_only_source_functional_obstruction"])
        self.assertFalse(self.acceptance["accepted_as_strict_source_law_or_coupling"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["role_bearing_ltotal_promoted"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["toe_closure_exported"])

    def test_written_documents_reference_guardrail(self):
        self.assertIn("P2810/S1760", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2810/S1760", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2810/S1760", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2810", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
