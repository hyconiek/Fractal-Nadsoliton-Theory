import json
import subprocess
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "fundamental_action_reconstruction"
SCRIPT = ROOT / "p2777_s1727_symmetry_source_selector_geometry_audit.py"
OUT = ROOT / "generated" / "p2777_s1727_symmetry_source_selector_geometry_audit.json"
MD = ROOT / "generated" / "p2777_s1727_symmetry_source_selector_geometry_audit.md"


class P2777SymmetrySourceSelectorGeometryAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(["python", str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2777_SYMMETRY_SOURCE_SELECTOR_GEOMETRY_AUDIT_NO_CLOSURE")
        self.assertEqual(self.payload["input_statuses"]["P2776"], "P2776_SMALL_GRAPH_FULL_SPECTRUM_UNIQUENESS_AUDIT_NO_CLOSURE")
        self.assertIn("symmetry", self.payload["audited_question"])

    def test_aut_z12_symmetry_does_not_select_orientation(self):
        selector = self.payload["symmetry_source_witness"]["selector_symmetry_witness"]
        self.assertEqual(selector["aut_units"], [1, 5, 7, 11])
        self.assertEqual(selector["orbit_of_directed_candidates"], [1, 5, 7, 11])
        self.assertTrue(selector["inversion_unit_present"])
        self.assertTrue(selector["plus_maps_to_minus_under_inversion"])
        self.assertTrue(selector["minus_maps_to_plus_under_inversion"])
        self.assertEqual(selector["aut_invariant_singleton_selector_count"], 0)
        self.assertFalse(selector["symmetry_selects_orientation"])

    def test_geometry_symmetry_is_pair_local_only(self):
        geometry = self.payload["symmetry_source_witness"]["geometry_symmetry_witness"]
        self.assertTrue(geometry["all_pair_geometries_vertex_transitive"])
        self.assertFalse(geometry["vertex_transitivity_selects_unique_geometry"])
        self.assertTrue(geometry["max_symmetry_selects_unique_geometry_on_pair"])
        self.assertEqual(geometry["max_automorphism_count_geometries"], ["torus_4x4"])
        rows = {row["geometry"]: row for row in geometry["geometry_rows"]}
        self.assertEqual(rows["torus_4x4"]["automorphism_count"], 384)
        self.assertEqual(rows["circulant_pm1_pm2"]["automorphism_count"], 32)
        self.assertTrue(rows["torus_4x4"]["vertex_transitive"])
        self.assertTrue(rows["circulant_pm1_pm2"]["vertex_transitive"])

    def test_acceptance_blocks_selector_and_geometry_closure(self):
        acceptance = self.payload["acceptance_matrix"]
        self.assertFalse(acceptance["accepted_as_selector_source_theorem"])
        self.assertTrue(acceptance["accepted_as_pair_local_max_symmetry_geometry_discriminator"])
        self.assertFalse(acceptance["accepted_as_canonical_nadsoliton_geometry_source"])
        self.assertFalse(acceptance["accepted_as_ltotal_or_toe_promotion"])
        self.assertIn("aut_z12_symmetry_selects_orientation", acceptance["missing_criteria"])
        self.assertIn("strict_source_law_for_max_symmetry_exported", acceptance["missing_criteria"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_documentation_and_recommendation(self):
        self.assertIn("symmetry-breaking/chiral/pseudoscalar source", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2697-P2777", self.payload["decision"]["next_honest_step"])
        self.assertIn("P2777/S1727", MD.read_text(encoding="utf-8"))
        self.assertIn("P2777/S1727", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2777/S1727", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2777", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
