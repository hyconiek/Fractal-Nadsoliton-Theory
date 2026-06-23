import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3048_s1998_memory_phase_area_odd_source_candidate.py"
OUT = ROOT / "generated" / "p3048_s1998_memory_phase_area_odd_source_candidate.json"
MD = ROOT / "generated" / "p3048_s1998_memory_phase_area_odd_source_candidate.md"

class P3048MemoryPhaseAreaOddSourceCandidateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3048_MEMORY_PHASE_AREA_ODD_SOURCE_CANDIDATE_BOUNDED_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3047"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["area_rows"], 6)
        self.assertEqual(cert["finite_nonzero_area_rows"], 5)
        self.assertEqual(cert["inversion_odd_rows"], 6)
        self.assertEqual(cert["accepted_strict_odd_source_value_rows"], 0)
        self.assertEqual(cert["source_acceptance_criteria"], 7)
        self.assertEqual(cert["satisfied_source_acceptance_criteria"], 3)
        self.assertFalse(cert["strict_inversion_odd_source_value_exported"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "MemoryPhaseSpaceArea_OddSourceCandidateMatrix")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["p3047_missing_odd_value_targeted"])
        self.assertTrue(obligations["three_point_area_not_pair_commutator_replay"])
        self.assertTrue(obligations["finite_nonzero_odd_value"])
        self.assertFalse(obligations["strict_nadsoliton_source_law"])
        self.assertFalse(obligations["nonconventional_orientation_and_coupling"])
        self.assertFalse(obligations["selector_readout_or_ltotal_installation"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3048/S1998", MD.read_text(encoding="utf-8"))
        self.assertIn("P3048/S1998", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3048/S1998", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3048", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
