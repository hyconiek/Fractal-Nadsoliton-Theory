import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3160_s2110_alpha_geo_phase_locking_no_root_of_unity_theorem.py"
OUT = ROOT / "generated" / "p3160_s2110_alpha_geo_phase_locking_no_root_of_unity_theorem.json"
MD = ROOT / "generated" / "p3160_s2110_alpha_geo_phase_locking_no_root_of_unity_theorem.md"


class P3160AlphaGeoPhaseLockingNoRootOfUnityTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_exact_theorem_shape(self):
        self.assertEqual(self.payload["status"], "P3160_ALPHA_GEO_PHASE_LOCKING_NO_ROOT_OF_UNITY_THEOREM")
        steps = self.payload["constructed_theoretical_objects"]["exact_theorem_steps"]
        self.assertEqual(len(steps), 5)
        self.assertIn("Gelfond-Schneider", steps[2]["claim"])
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["exact_irrationality_theorem_exported"])
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["no_finite_ZN_phase_locking_proved"])

    def test_computation_rows_and_no_exports(self):
        self.assertEqual(len(self.payload["denominator_scan_rows_N_1_to_144"]), 144)
        self.assertEqual(len(self.payload["continued_fraction_ladder"]), 11)
        self.assertEqual(self.payload["finite_theorem"]["constructed_objects"], 5)
        self.assertEqual(self.payload["finite_theorem"]["accepted_closure_exports"], 0)
        self.assertTrue(all(row["exact_locking_excluded_by_theorem"] for row in self.payload["denominator_scan_rows_N_1_to_144"]))
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))
        self.assertIn("Do not continue alpha_geo/pi phase-locking", self.payload["decision"]["next_honest_step"])

    def test_docs_updated(self):
        self.assertIn("P3160/S2110", MD.read_text(encoding="utf-8"))
        self.assertIn("P3160/S2110", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3160/S2110", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3160", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
