import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2778_s1728_max_symmetry_16node_geometry_source_audit.py"
OUT = ROOT / "generated" / "p2778_s1728_max_symmetry_16node_geometry_source_audit.json"
MD = ROOT / "generated" / "p2778_s1728_max_symmetry_16node_geometry_source_audit.md"


class P2778MaxSymmetry16NodeGeometrySourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2778_MAX_SYMMETRY_16NODE_GEOMETRY_SOURCE_AUDIT_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2777"], "P2777_SYMMETRY_SOURCE_SELECTOR_GEOMETRY_AUDIT_NO_CLOSURE")
        self.assertIn("maximal-automorphism", self.payload["audited_question"])

    def test_max_symmetry_fails_to_select_torus_on_declared_class(self):
        witness = self.payload["max_symmetry_16node_witness"]
        self.assertEqual(witness["candidate_count"], 19)
        self.assertTrue(witness["all_candidates_connected_4_regular"])
        self.assertTrue(witness["all_candidates_vertex_transitive"])
        self.assertEqual(witness["max_automorphism_count"], 4096)
        self.assertEqual(witness["max_geometry_count"], 2)
        self.assertEqual(witness["max_geometry_names"], ["circulant_pm1_pm7", "circulant_pm3_pm5"])
        self.assertEqual(witness["torus_4x4_automorphism_count"], 384)
        self.assertFalse(witness["max_symmetry_selects_torus_4x4"])
        self.assertFalse(witness["max_symmetry_selects_unique_labeled_geometry"])

    def test_acceptance_blocks_max_symmetry_geometry_source(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertTrue(acceptance["facts"]["strict_candidate_law_supplied"])
        self.assertTrue(acceptance["facts"]["declared_16node_class_audited"])
        self.assertFalse(acceptance["facts"]["max_symmetry_selects_torus_4x4"])
        self.assertFalse(acceptance["facts"]["max_symmetry_selects_unique_labeled_geometry"])
        self.assertFalse(acceptance["accepted_as_max_symmetry_geometry_source"])
        self.assertFalse(acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("max_symmetry_selects_torus_4x4", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("different sourced symmetry functional", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2778", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2778/S1728", MD.read_text(encoding="utf-8"))
        self.assertIn("P2778/S1728", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2778/S1728", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2778", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
