import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3138_s2088_nonfourier_dhl_extremum_joint_source_audit.py"
OUT = ROOT / "generated" / "p3138_s2088_nonfourier_dhl_extremum_joint_source_audit.json"
MD = ROOT / "generated" / "p3138_s2088_nonfourier_dhl_extremum_joint_source_audit.md"


class P3138NonFourierDHLExtremumJointSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_constructed_object(self):
        self.assertEqual(self.payload["status"], "P3138_NONFOURIER_E_DHL_EXTREMUM_JOINT_SOURCE_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3137"])
        self.assertIn("No Fourier coefficients", self.payload["constructed_object"]["why_non_fourier"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["profiles_tested"], 120)
        self.assertEqual(cert["translation_t1_covariant_receiver_rows"], 120)
        self.assertGreater(cert["inversion_paired_rows"], 0)
        self.assertGreater(cert["multi_tie_positive_slope_rows"], 0)
        self.assertEqual(cert["accepted_import_free_joint_sources"], 0)
        self.assertLess(cert["gates_passed"], cert["gates_required"])

    def test_gate_rows_and_exports(self):
        gates = {row["gate"]: row["passed"] for row in self.payload["source_gate_rows"]}
        self.assertTrue(gates["non_fourier_formula"])
        self.assertTrue(gates["computes_nonempty_origin_receiver"])
        self.assertFalse(gates["translation_quotient_invariant"])
        self.assertFalse(gates["inversion_polarity_unpaired"])
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))

    def test_docs_updated(self):
        self.assertIn("P3138/S2088", MD.read_text(encoding="utf-8"))
        self.assertIn("P3138/S2088", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3138/S2088", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3138", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
