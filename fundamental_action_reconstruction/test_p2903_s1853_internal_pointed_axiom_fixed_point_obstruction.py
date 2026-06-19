import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2903_s1853_internal_pointed_axiom_fixed_point_obstruction.py"
JSON_PATH = ROOT / "generated" / "p2903_s1853_internal_pointed_axiom_fixed_point_obstruction.json"
MD_PATH = ROOT / "generated" / "p2903_s1853_internal_pointed_axiom_fixed_point_obstruction.md"


class P2903InternalPointedAxiomFixedPointObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.objects = cls.payload["constructed_theoretical_objects"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(cls_status := self.payload["status"], "P2903_INTERNAL_POINTED_AXIOM_FIXED_POINT_OBSTRUCTION_NO_STRICT_SOURCE")
        self.assertIn("P2903", cls_status)
        self.assertTrue(self.acceptance["p2902_rechecked"])

    def test_fixed_point_obstruction_counts(self):
        self.assertEqual(self.acceptance["target_point_count"], 24)
        self.assertEqual(self.acceptance["translation_orbit_count"], 2)
        self.assertEqual(self.acceptance["translation_orbit_sizes"], [12, 12])
        self.assertEqual(self.acceptance["fixed_point_count"], 0)
        self.assertEqual(self.acceptance["equivariant_map_count_from_translation_trivial_source"], 0)
        self.assertFalse(self.acceptance["internal_translation_neutral_source_can_select_pointed_axiom"])
        self.assertTrue(self.acceptance["requires_translation_breaking_strict_source"])

    def test_constructed_proof_objects(self):
        self.assertEqual(self.objects["finite_target_point_count"], 24)
        self.assertEqual(self.objects["fixed_points"], [])
        self.assertEqual(self.objects["equivariant_maps_from_trivial_source"], [])
        self.assertEqual(len(self.objects["proof_steps"]), 5)

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.flags["internal_pointed_signed_axiom_exported"])
        self.assertFalse(self.flags["nonproxy_ltotal_exported"])
        self.assertFalse(self.flags["toe_closure_exported"])

    def test_documents_updated(self):
        self.assertIn("P2903/S1853", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2903/S1853", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2903/S1853", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2903", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
