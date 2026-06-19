import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2933_s1883_aut_breaking_prime_coordinate_source_space_decomposition.py"
OUT = ROOT / "generated" / "p2933_s1883_aut_breaking_prime_coordinate_source_space_decomposition.json"
MD = ROOT / "generated" / "p2933_s1883_aut_breaking_prime_coordinate_source_space_decomposition.md"


class P2933AutBreakingPrimeCoordinateSourceSpaceDecompositionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2933_AUT_BREAKING_PRIME_COORDINATE_SOURCE_SPACE_DECOMPOSITION_NO_SOURCE_LAW")
        self.assertIsNotNone(self.payload["input_hashes"]["P2932"])

    def test_decomposition_certificate(self):
        cert = self.payload["decomposition_certificate"]
        self.assertEqual(cert["product_equation_count"], 30)
        self.assertEqual(cert["aut_equation_count"], 44)
        self.assertEqual(cert["product_rank"], 6)
        self.assertEqual(cert["product_nullity"], 5)
        self.assertEqual(cert["combined_rank"], 11)
        self.assertEqual(cert["combined_nullity"], 0)
        self.assertEqual(cert["prime_coordinate_basis_count"], 5)
        self.assertTrue(cert["all_basis_vectors_break_aut_invariance"])
        self.assertEqual(cert["accepted_candidate_count"], 0)

    def test_future_object_and_no_closure_flags(self):
        source_space = self.payload["constructed_theoretical_objects"]["source_space"]
        self.assertEqual(source_space["name"], "AutBreaking_PrimeCoordinate_Source_Space")
        self.assertEqual(source_space["aut_breaking_quotient_dimension"], 5)
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_prime_log_value_source_exported", "strict_delta_eta_source_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2933/S1883", MD.read_text(encoding="utf-8"))
        self.assertIn("P2933/S1883", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2933/S1883", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2933", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
