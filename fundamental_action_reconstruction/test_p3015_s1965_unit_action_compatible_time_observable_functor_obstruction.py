import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3015_s1965_unit_action_compatible_time_observable_functor_obstruction.py"
OUT = ROOT / "generated" / "p3015_s1965_unit_action_compatible_time_observable_functor_obstruction.json"
MD = ROOT / "generated" / "p3015_s1965_unit_action_compatible_time_observable_functor_obstruction.md"

class P3015UnitActionCompatibleTimeObservableFunctorTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3015_UNIT_ACTION_COMPATIBLE_TIME_OBSERVABLE_FUNCTOR_OBSTRUCTION_NO_CLOSURE")
        self.assertIsNotNone(self.payload["input_hashes"]["P3014"])

    def test_unit_compatibility_and_successor_obstruction(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["orbit_count"], 6)
        self.assertEqual(cert["unit_compatibility_rows"], 48)
        self.assertEqual(cert["unit_compatibility_failure_count"], 0)
        self.assertEqual(cert["successor_row_count"], 6)
        self.assertGreater(cert["bad_successor_orbit_count"], 0)
        self.assertFalse(cert["accepted_as_time_observable_generator"])

    def test_functor_obligations_and_exports(self):
        functor = self.payload["constructed_theoretical_objects"]["functor"]
        self.assertEqual(functor["candidate"], "StrictKernelOrbitQuotient_TimeObservableFunctor")
        obligations = {row["obligation"]: row["satisfied"] for row in functor["proof_obligations"]}
        self.assertTrue(obligations["explicit_input_output_types"])
        self.assertTrue(obligations["observer_independent_formula"])
        self.assertTrue(obligations["unit_action_compatible"])
        self.assertFalse(obligations["well_defined_clock_successor"])
        self.assertFalse(obligations["directed_time_arrow_without_selector_import"])
        self.assertFalse(obligations["eom_hamiltonian_installation"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3015/S1965", MD.read_text(encoding="utf-8"))
        self.assertIn("P3015/S1965", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3015/S1965", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3015", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
