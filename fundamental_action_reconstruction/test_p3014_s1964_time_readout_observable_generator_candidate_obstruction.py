import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3014_s1964_time_readout_observable_generator_candidate_obstruction.py"
OUT = ROOT / "generated" / "p3014_s1964_time_readout_observable_generator_candidate_obstruction.json"
MD = ROOT / "generated" / "p3014_s1964_time_readout_observable_generator_candidate_obstruction.md"

class P3014TimeReadoutObservableGeneratorTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3014_TIME_READOUT_OBSERVABLE_GENERATOR_CANDIDATE_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3013"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertTrue(cert["explicit_io"])
        self.assertTrue(cert["observer_independent_formula"])
        self.assertEqual(cert["unit_equivariance_rows"], 48)
        self.assertGreater(cert["unit_equivariance_failure_count"], 0)
        self.assertFalse(cert["accepted_as_time_observable_generator"])

    def test_witness_obligations_and_exports(self):
        witness = self.payload["constructed_theoretical_objects"]["witness"]
        self.assertEqual(witness["candidate"], "StrictKernelFiniteDifference_TimeObservableGenerator")
        self.assertFalse(witness["accepted_as_time_observable_generator"])
        obligations = {row["obligation"]: row["satisfied"] for row in witness["proof_obligations"]}
        self.assertTrue(obligations["explicit_input_output_types"])
        self.assertTrue(obligations["observer_independent_formula"])
        self.assertFalse(obligations["unit_action_compatible"])
        self.assertFalse(obligations["eom_hamiltonian_installation"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3014/S1964", MD.read_text(encoding="utf-8"))
        self.assertIn("P3014/S1964", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3014/S1964", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3014", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
