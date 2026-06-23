import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3059_s2009_polarity_odd_source_law_synthesis_verifier.py"
OUT = ROOT / "generated" / "p3059_s2009_polarity_odd_source_law_synthesis_verifier.json"
MD = ROOT / "generated" / "p3059_s2009_polarity_odd_source_law_synthesis_verifier.md"

class P3059PolarityOddSourceLawSynthesisVerifierTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3059_POLARITY_ODD_SOURCE_LAW_SYNTHESIS_SIGN_PAIR_OBSTRUCTION_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3058"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["odd_basis_size"], 3)
        self.assertEqual(cert["coefficient_values_per_basis"], 5)
        self.assertEqual(cert["total_nonzero_coefficient_vectors"], 124)
        self.assertEqual(cert["nonzero_signed_candidates"], 106)
        self.assertEqual(cert["zero_value_candidates"], 18)
        self.assertEqual(cert["positive_polarity_candidates"], 53)
        self.assertEqual(cert["negative_polarity_candidates"], 53)
        self.assertEqual(cert["sign_pair_orbits"], 53)
        self.assertEqual(cert["accepted_unique_boundary_conditions"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 3)

    def test_synthesis_module_and_pairs(self):
        obj = self.payload["constructed_theoretical_objects"]
        module = obj["synthesis_module"]
        self.assertEqual(module["object"], "StrictPolarityOddSourceLawBoundaryConditionSynthesisModule")
        self.assertEqual(module["carrier_symbol"], "G_selector")
        self.assertEqual(module["target_atom"], "strict_polarity_odd_source_law_boundary_condition")
        self.assertEqual(len(module["basis_clues"]), 3)
        for row in obj["candidate_rows"]:
            self.assertTrue(row["inversion_odd"])
            self.assertFalse(row["accepted_as_unique_nonpremise_boundary_condition"])
            self.assertEqual(row["paired_base_value"], -row["base_value"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3059/S2009", MD.read_text(encoding="utf-8"))
        self.assertIn("P3059/S2009", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3059/S2009", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3059", (REPO / "AGENTS.md").read_text(encoding="utf-8"))
        self.assertIn("`106`", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("`53` positive", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
