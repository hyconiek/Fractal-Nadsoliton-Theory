import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2906_s1856_xi_strict_asymmetry_theorem_acceptance_obstruction.py"
JSON_PATH = ROOT / "generated" / "p2906_s1856_xi_strict_asymmetry_theorem_acceptance_obstruction.json"
MD_PATH = ROOT / "generated" / "p2906_s1856_xi_strict_asymmetry_theorem_acceptance_obstruction.md"


class P2906XiStrictAsymmetryTheoremAcceptanceObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
        cls.payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
        cls.acceptance = cls.payload["acceptance_matrix"]
        cls.objects = cls.payload["constructed_theoretical_objects"]
        cls.flags = cls.payload["decision"]["negative_export_flags"]

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2906_XI_STRICT_ASYMMETRY_THEOREM_ACCEPTANCE_OBSTRUCTION_NO_EXPORT")
        self.assertTrue(self.acceptance["p2905_rechecked_no_positive_provenance"])

    def test_target_orbits_and_fixed_points(self):
        self.assertEqual(self.acceptance["xi_target_count"], 24)
        self.assertEqual(self.acceptance["translation_orbit_count"], 2)
        self.assertEqual(self.acceptance["translation_orbit_sizes"], [12, 12])
        self.assertEqual(self.acceptance["translation_fixed_point_count"], 0)

    def test_sign_and_origin_obligations_separated(self):
        self.assertEqual(self.acceptance["chiral_sign_only_remaining_basepoints"], 12)
        self.assertEqual(self.acceptance["origin_only_remaining_polarities"], 2)
        self.assertTrue(self.acceptance["joint_origin_sign_theorem_required"])
        self.assertFalse(self.acceptance["joint_origin_sign_theorem_exported"])
        rows = self.objects["acceptance_obstruction_rows"]
        self.assertEqual(len(rows), 4)
        self.assertIn("chiral/sign source only", rows[1]["candidate_source_kind"])

    def test_false_export_flags(self):
        self.assertFalse(any(self.flags.values()))
        self.assertFalse(self.flags["xi_0_plus_strict_source_exported"])
        self.assertFalse(self.flags["nonproxy_ltotal_exported"])
        self.assertFalse(self.flags["toe_closure_exported"])

    def test_documents_updated(self):
        self.assertIn("P2906/S1856", MD_PATH.read_text(encoding="utf-8"))
        self.assertIn("P2906/S1856", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2906/S1856", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2906", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
