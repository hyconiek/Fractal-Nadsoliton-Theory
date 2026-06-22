import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3020_s1970_clock_unit_theorem_candidate_obstruction.py"
OUT = ROOT / "generated" / "p3020_s1970_clock_unit_theorem_candidate_obstruction.json"
MD = ROOT / "generated" / "p3020_s1970_clock_unit_theorem_candidate_obstruction.md"

class P3020ClockUnitTheoremCandidateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3020_CLOCK_UNIT_THEOREM_CANDIDATE_OBSTRUCTION_NO_CLOSURE")
        self.assertIsNotNone(self.payload["input_hashes"]["P3019"])

    def test_finite_clock_candidates_but_no_clock_unit(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["clock_candidate_count"], 3)
        self.assertEqual(cert["positive_clock_candidates"], 3)
        self.assertEqual(cert["unit_count"], 4)
        self.assertEqual(cert["plus_one_tick_preserving_unit_count"], 1)
        self.assertTrue(cert["dominant_dft_modes"])
        self.assertFalse(cert["accepted_as_strict_clock_unit_theorem"])

    def test_obligations_and_negative_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "ClockUnitTheoremCandidate_UnitActionSpectralObstructionMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["typed_time_observable_input"])
        self.assertTrue(obligations["positive_finite_clock_candidates"])
        self.assertFalse(obligations["unit_invariant_directed_tick"])
        self.assertFalse(obligations["local_successor_from_cycle_period"])
        self.assertFalse(obligations["physical_frequency_or_clock_unit_source"])
        self.assertFalse(obligations["action_hamiltonian_coupling"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3020/S1970", MD.read_text(encoding="utf-8"))
        self.assertIn("P3020/S1970", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3020/S1970", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3020", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
