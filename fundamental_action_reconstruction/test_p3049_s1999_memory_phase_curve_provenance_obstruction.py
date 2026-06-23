import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3049_s1999_memory_phase_curve_provenance_obstruction.py"
OUT = ROOT / "generated" / "p3049_s1999_memory_phase_curve_provenance_obstruction.json"
MD = ROOT / "generated" / "p3049_s1999_memory_phase_curve_provenance_obstruction.md"

class P3049MemoryPhaseCurveProvenanceObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3049_MEMORY_PHASE_CURVE_PROVENANCE_OBSTRUCTION_BOUNDED_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3048"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["provenance_rows"], 5)
        self.assertEqual(cert["strict_source_provenance_rows"], 0)
        self.assertEqual(cert["aut_area_rows"], 24)
        self.assertEqual(cert["aut_equivariance_verified_rows"], 24)
        self.assertEqual(cert["orientation_reversing_rows"], 12)
        self.assertEqual(cert["source_acceptance_criteria"], 8)
        self.assertEqual(cert["satisfied_source_acceptance_criteria"], 3)
        self.assertFalse(cert["strict_KM_phase_curve_source_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "MemoryPhaseCurve_ProvenanceObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["p3048_phase_curve_provenance_targeted"])
        self.assertTrue(obligations["source_channels_explicitly_separated"])
        self.assertTrue(obligations["finite_aut_relabeling_test"])
        self.assertFalse(obligations["strict_KM_source_provenance"])
        self.assertFalse(obligations["chart_independent_orientation"])
        self.assertFalse(obligations["selector_readout_or_ltotal_installation"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3049/S1999", MD.read_text(encoding="utf-8"))
        self.assertIn("P3049/S1999", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3049/S1999", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3049", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
