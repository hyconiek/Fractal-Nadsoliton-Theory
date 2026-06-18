import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2846_s1796_label_safe_vertex_localization_pullback_candidate_audit.py"
JSON_PATH = ROOT / "generated" / "p2846_s1796_label_safe_vertex_localization_pullback_candidate_audit.json"
MD_PATH = ROOT / "generated" / "p2846_s1796_label_safe_vertex_localization_pullback_candidate_audit.md"


class P2846LabelSafeVertexLocalizationPullbackAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(
            self.payload["status"],
            "P2846_LABEL_SAFE_VERTEX_LOCALIZATION_PULLBACK_CANDIDATE_NO_GO_NO_CLOSURE",
        )
        audit = self.payload["label_safe_vertex_localization_pullback_audit"]
        self.assertEqual(
            audit["input_statuses_rechecked"]["P2845"],
            "P2845_UNIT_BEARING_TYPED_SOURCE_COUPLING_DIMENSIONAL_OBSTRUCTION_NO_CLOSURE",
        )

    def test_full_carrier_and_profiles(self):
        audit = self.payload["label_safe_vertex_localization_pullback_audit"]
        self.assertEqual(audit["carrier_check"]["decoded_graph_count"], 16828)
        self.assertTrue(audit["carrier_check"]["coverage_ok"])
        rows = audit["candidate_rows"]
        self.assertEqual(len(rows), 5)
        self.assertFalse(rows["uniform_vertex_measure"]["profile_stats"]["nonconstant_on_full_carrier"])
        self.assertFalse(rows["triangle_count_vertex_measure"]["profile_stats"]["nonconstant_on_full_carrier"])
        self.assertGreater(rows["four_cycle_count_vertex_measure"]["profile_stats"]["distinct_profile_count"], 1)
        self.assertGreater(rows["local_motif_wl_vertex_measure"]["profile_stats"]["distinct_profile_count"], 1)

    def test_all_candidates_rejected_as_strict_pullback(self):
        audit = self.payload["label_safe_vertex_localization_pullback_audit"]
        self.assertGreater(audit["nonconstant_candidate_count"], 0)
        self.assertEqual(audit["accepted_candidate_count"], 0)
        for row in audit["candidate_rows"].values():
            self.assertTrue(row["premises"]["label_gauge_invariant"])
            self.assertTrue(row["premises"]["normalized_finite_vertex_measure"])
            self.assertFalse(row["premises"]["canonical_vertex_to_field_support"])
            self.assertFalse(row["premises"]["spacetime_pullback_formula"])
            self.assertFalse(row["accepted_as_strict_localization_pullback"])

    def test_acceptance_and_negative_exports(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["accepted_as_label_safe_localization_obstruction_audit"])
        self.assertFalse(acceptance["exports_strict_localization_pullback"])
        flags = self.payload["decision"]["negative_export_flags"]
        self.assertFalse(flags["strict_localization_pullback_exported"])
        self.assertFalse(flags["nonproxy_ltotal_term_exported"])
        self.assertFalse(flags["eom_closure_exported"])
        self.assertFalse(flags["hamiltonian_closure_exported"])
        self.assertFalse(flags["toe_closure_exported"])

    def test_documents_updated(self):
        self.assertIn("P2846/S1796", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2846/S1796", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2846/S1796", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2846", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
