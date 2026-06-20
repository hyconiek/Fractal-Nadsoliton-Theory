import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2962_s1912_kc_exchange_symmetry_source_obstruction.py"
OUT = ROOT / "generated" / "p2962_s1912_kc_exchange_symmetry_source_obstruction.json"
MD = ROOT / "generated" / "p2962_s1912_kc_exchange_symmetry_source_obstruction.md"

class P2962KCExchangeSymmetrySourceObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2962_KC_EXCHANGE_SYMMETRY_SOURCE_OBSTRUCTION_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2961"])

    def test_artifact_invariants_block_exchange(self):
        cert = self.payload["obstruction_certificate"]
        self.assertTrue(cert["coefficient_quotient_candidate_available"])
        self.assertFalse(cert["support_cardinality_match"])
        self.assertFalse(cert["nonzero_multiset_match"])
        self.assertFalse(cert["typed_provenance_equivalence_exported"])
        self.assertFalse(cert["strict_KC_exchange_source_exported"])
        invariants = {r["invariant"]: r for r in self.payload["constructed_theoretical_objects"]["artifact_invariant_rows"]}
        self.assertEqual(invariants["support_cardinality"]["K"], 2)
        self.assertEqual(invariants["support_cardinality"]["C"], 3)
        self.assertEqual(invariants["coordinate_permutation_to_counterpart"]["K_to_C_maps"], 0)

    def test_obligations_and_acceptance_matrix(self):
        obligations = {r["obligation"]: r["satisfied"] for r in self.payload["constructed_theoretical_objects"]["proof_obligation_rows"]}
        self.assertTrue(obligations["P2961_coefficient_quotient_exchange_candidate_available"])
        self.assertFalse(obligations["K_and_C_support_cardinalities_match"])
        self.assertFalse(obligations["K_and_C_nonzero_multisets_match"])
        self.assertFalse(obligations["strict_nadsoliton_artifact_symmetry_source_exported"])
        matrix = self.payload["constructed_theoretical_objects"]["finite_acceptance_matrix"]
        self.assertEqual(len(matrix), 64)
        self.assertEqual(sum(1 for row in matrix if row["accepts_strict_KC_exchange_source"]), 1)

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2962/S1912", MD.read_text(encoding="utf-8"))
        self.assertIn("P2962/S1912", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2962/S1912", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2962", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
