import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2991_s1941_z12_additive_fourier_character_source_candidate_obstruction.py"
OUT = ROOT / "generated" / "p2991_s1941_z12_additive_fourier_character_source_candidate_obstruction.json"
MD = ROOT / "generated" / "p2991_s1941_z12_additive_fourier_character_source_candidate_obstruction.md"

class P2991Z12AdditiveFourierCharacterSourceCandidateObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2991_Z12_ADDITIVE_FOURIER_CHARACTER_SOURCE_CANDIDATE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P2990"])

    def test_character_certificate(self):
        cert = self.payload["character_certificate"]
        self.assertEqual(cert["character_count"], 12)
        self.assertEqual(cert["orthogonality_matrix_rows"], 144)
        self.assertTrue(cert["all_orthogonality_checks_pass"])
        self.assertEqual(cert["conductor_histogram"], {"1": 1, "2": 1, "3": 2, "4": 2, "6": 2, "12": 4})
        self.assertEqual(cert["accepted_strict_character_sources"], [])
        self.assertEqual(cert["acceptance_matrix_rows"], 256)
        self.assertEqual(cert["accepted_rows"], 1)

    def test_obligations_and_rows(self):
        obj = self.payload["constructed_theoretical_objects"]
        obligations = {r["obligation"]: r["satisfied"] for r in obj["proof_obligation_rows"]}
        self.assertTrue(obligations["finite_character_table_enumerated"])
        self.assertTrue(obligations["full_orthogonality_matrix_verified"])
        self.assertTrue(obligations["conductor_kernel_classification"])
        self.assertFalse(obligations["nonpremise_frequency_selector"])
        self.assertFalse(obligations["strict_nadsoliton_character_provenance"])
        self.assertFalse(obligations["source_coupling_theorem"])
        self.assertFalse(obligations["unit_bearing_measure_or_action_density"])
        self.assertFalse(obligations["accepted_current_strict_character_source"])
        witness = obj["fourier_character_witness"]
        self.assertEqual(len(witness["character_rows"]), 12)
        self.assertEqual(len(witness["orthogonality_matrix_rows"]), 144)
        self.assertFalse(any(r["accepted_strict_character_source"] for r in witness["character_rows"]))

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2991/S1941", MD.read_text(encoding="utf-8"))
        self.assertIn("P2991/S1941", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2991/S1941", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2991", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
