import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3025_s1975_dissipation_physical_tick_hamiltonian_unit_obstruction.py"
OUT = ROOT / "generated" / "p3025_s1975_dissipation_physical_tick_hamiltonian_unit_obstruction.json"
MD = ROOT / "generated" / "p3025_s1975_dissipation_physical_tick_hamiltonian_unit_obstruction.md"

class P3025DissipationPhysicalTickHamiltonianUnitTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3025_DISSIPATION_PHYSICAL_TICK_HAMILTONIAN_UNIT_OBSTRUCTION_NO_CLOSURE")
        self.assertIsNotNone(self.payload["input_hashes"]["P3024"])

    def test_tick_candidates_positive_but_not_physical_units(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["tick_candidate_count"], 4)
        self.assertEqual(cert["positive_tick_candidates"], 4)
        self.assertEqual(cert["rescaling_independent_tick_rows"], 1)
        self.assertEqual(cert["hamiltonian_unit_rows"], 0)
        self.assertFalse(cert["accepted_as_physical_tick_hamiltonian_unit_theorem"])

    def test_obligations_and_rescaling_rows(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "DissipationChainPhysicalTickHamiltonianUnit_RescalingCouplingObstructionMatrix")
        tick_rescalings = {row["candidate_tick"]: row["tick_rescaling"] for row in obj["tick_rows"]}
        self.assertEqual(tick_rescalings["label_tick"], "c^0")
        self.assertEqual(tick_rescalings["mean_positive_drop_tick"], "c^1")
        self.assertEqual(tick_rescalings["rms_drop_tick"], "c^1")
        self.assertEqual(tick_rescalings["inverse_mean_slope_tick"], "c^-1")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["attacks_single_P3024_remaining_theorem"])
        self.assertTrue(obligations["positive_tick_candidates_constructed"])
        self.assertFalse(obligations["observable_rescaling_independent_tick"])
        self.assertFalse(obligations["independent_action_quantum_source"])
        self.assertFalse(obligations["physical_Hamiltonian_unit_coupling"])
        self.assertFalse(obligations["strict_physical_tick_theorem_exported"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3025/S1975", MD.read_text(encoding="utf-8"))
        self.assertIn("P3025/S1975", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3025/S1975", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3025", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
