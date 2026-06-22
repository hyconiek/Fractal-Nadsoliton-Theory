import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3019_s1969_observable_unit_readout_source_obstruction.py"
OUT = ROOT / "generated" / "p3019_s1969_observable_unit_readout_source_obstruction.json"
MD = ROOT / "generated" / "p3019_s1969_observable_unit_readout_source_obstruction.md"

class P3019ObservableUnitReadoutSourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3019_OBSERVABLE_UNIT_READOUT_SOURCE_OBSTRUCTION_NO_CLOSURE")
        self.assertIsNotNone(self.payload["input_hashes"]["P3018"])

    def test_positive_units_but_no_absolute_unit(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["candidate_unit_count"], 4)
        self.assertEqual(cert["positive_candidate_units"], 4)
        self.assertEqual(cert["scale_row_count"], 4)
        self.assertEqual(cert["all_units_rescale_rows"], 4)
        self.assertEqual(cert["absolute_unit_fixed_rows"], 1)
        self.assertFalse(cert["accepted_as_strict_observable_unit_source"])

    def test_obligations_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "ObservableUnitReadoutSource_RescalingOrbitObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["typed_time_observable_input"])
        self.assertTrue(obligations["positive_finite_candidate_units"])
        self.assertTrue(obligations["observer_independent_formula"])
        self.assertFalse(obligations["observable_rescaling_invariant_absolute_unit"])
        self.assertFalse(obligations["strict_physical_readout_unit_source"])
        self.assertFalse(obligations["clock_hamiltonian_coupling"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3019/S1969", MD.read_text(encoding="utf-8"))
        self.assertIn("P3019/S1969", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3019/S1969", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3019", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
