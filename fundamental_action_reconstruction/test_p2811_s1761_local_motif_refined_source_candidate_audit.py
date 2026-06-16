import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
JSON_PATH = ROOT / "generated" / "p2811_s1761_local_motif_refined_source_candidate_audit.json"
MD_PATH = ROOT / "generated" / "p2811_s1761_local_motif_refined_source_candidate_audit.md"


class P2811LocalMotifRefinedSourceCandidateAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(ROOT / "p2811_s1761_local_motif_refined_source_candidate_audit.py")], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.audit = cls.payload["local_motif_refined_audit"]
        cls.acceptance = cls.payload["acceptance_matrix"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2811_LOCAL_MOTIF_REFINED_SOURCE_CANDIDATE_AUDIT_REFINES_BUT_STILL_OBSTRUCTED_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2810"], "P2810_SPECTRAL_PAIR_ONLY_SOURCE_FUNCTIONAL_OBSTRUCTION_NO_STRICT_SOURCE_LAW_NO_CLOSURE")

    def test_exact_refinement_counts(self):
        self.assertEqual(self.audit["decoded_graph_count"], 16828)
        self.assertEqual(self.audit["spectral_pair_class_count_from_p2810"], 16211)
        self.assertEqual(self.audit["refined_class_count"], 16691)
        self.assertEqual(self.audit["canonical_target_class_count_from_p2810"], 16828)
        self.assertEqual(self.audit["refined_collision_class_count"], 132)
        self.assertEqual(self.audit["refined_collision_graph_count"], 269)
        self.assertEqual(self.audit["refined_max_class_size"], 3)
        self.assertEqual(self.audit["remaining_defect_canonical_minus_refined"], 137)
        self.assertEqual(self.audit["defect_reduction_vs_spectral_pair_only"], 480)

    def test_acceptance_boundaries(self):
        self.assertTrue(self.acceptance["accepted_as_local_motif_refinement_audit"])
        self.assertFalse(self.acceptance["accepted_as_complete_canonical_source_carrier"])
        self.assertFalse(self.acceptance["accepted_as_strict_source_law_or_coupling"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["role_bearing_ltotal_promoted"])
        self.assertFalse(self.payload["decision"]["negative_export_flags"]["toe_closure_exported"])

    def test_written_documents_reference_guardrail(self):
        self.assertIn("P2811/S1761", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2811/S1761", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2811/S1761", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2811", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
