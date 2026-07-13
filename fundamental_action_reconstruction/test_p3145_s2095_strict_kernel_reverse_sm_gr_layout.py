import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3145_s2095_strict_kernel_reverse_sm_gr_layout.py"
OUT = ROOT / "generated" / "p3145_s2095_strict_kernel_reverse_sm_gr_layout.json"
MD = ROOT / "generated" / "p3145_s2095_strict_kernel_reverse_sm_gr_layout.md"


class P3145StrictKernelReverseSmGrLayoutTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_counts(self):
        self.assertEqual(self.payload["status"], "P3145_STRICT_KERNEL_REVERSE_SM_GR_LAYOUT_MATRIX_RECEIVER_ONLY_NO_CLOSURE")
        self.assertTrue(self.payload["input_hashes"]["P3144"])
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["physical_properties_audited"], 10)
        self.assertEqual(counts["strict_carrier_rows"], 10)
        self.assertEqual(counts["receiver_scaffold_rows"], 10)
        self.assertEqual(counts["source_law_rows"], 0)
        self.assertEqual(counts["unit_or_calibration_rows"], 0)
        self.assertEqual(counts["nonproxy_eom_rows"], 0)
        self.assertEqual(counts["closed_rows"], 0)
        self.assertEqual(counts["receiver_layout_only_rows"], 10)

    def test_layout_rows_cover_sm_gr_properties(self):
        props = {row["physical_property"] for row in self.payload["layout_rows"]}
        self.assertIn("GR stress-energy and Einstein dynamics", props)
        self.assertIn("SM gauge scaffold SU(3)xSU(2)xU(1)", props)
        self.assertIn("fermion masses, Yukawa-like hierarchy, flavor", props)
        self.assertIn("dimensionful units: length, action, time, calibration", props)
        self.assertTrue(all(row["layout_status"] == "receiver_layout_only" for row in self.payload["layout_rows"]))

    def test_decision_preserves_no_closure(self):
        decision = self.payload["decision"]
        self.assertIn("receiver architecture", decision["bounded_result"])
        self.assertIn("Continue axiomatizing only if", decision["axiomatic_route_recommendation"])
        self.assertIn("do not claim SM/GR reduction or ToE", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_updated(self):
        self.assertIn("P3145/S2095", MD.read_text(encoding="utf-8"))
        self.assertIn("P3145/S2095", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3145/S2095", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3145", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
