import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2958_s1908_p2938_strict_provenance_theorem_interface.py"
OUT = ROOT / "generated" / "p2958_s1908_p2938_strict_provenance_theorem_interface.json"
MD = ROOT / "generated" / "p2958_s1908_p2938_strict_provenance_theorem_interface.md"


class P2958P2938StrictProvenanceTheoremInterfaceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2958_P2938_STRICT_PROVENANCE_THEOREM_INTERFACE_NO_STRICT_EXPORT")
        for key in ["P2938", "P2953", "P2957"]:
            self.assertIsNotNone(self.payload["input_hashes"][key])

    def test_decomposition_square(self):
        square = self.payload["constructed_theoretical_objects"]["finite_provenance_square"]
        self.assertEqual(square["kernel_excess_vector_K"], [1, 2, 0, 0, 0])
        self.assertEqual(square["character_negativity_vector_C"], [0, 0, 2, 2, 2])
        self.assertEqual(square["target_vector_V"], [1, 2, 2, 2, 2])
        self.assertTrue(square["K_plus_C_equals_V"])
        self.assertEqual(square["sum_V"], 9)
        self.assertEqual(square["eta_from_sum_over_five_primes"]["as_string"], "9/5")
        self.assertEqual(square["support_overlap_count"], 0)

    def test_provenance_obligations_and_matrix(self):
        obligations = {row["obligation"]: row for row in self.payload["constructed_theoretical_objects"]["provenance_obligation_rows"]}
        self.assertTrue(obligations["finite_U12_character_aggregate_constructed"]["satisfied"])
        self.assertTrue(obligations["product_additive_carrier_available"]["satisfied"])
        self.assertFalse(obligations["strict_nadsoliton_functor_to_U12_aggregate_exported"]["satisfied"])
        self.assertFalse(obligations["nonconventional_equal_summand_law_exported"]["satisfied"])
        self.assertFalse(obligations["aggregate_localizer_not_just_finite_carrier"]["satisfied"])
        self.assertFalse(obligations["ratio_package_beta_unit_coupling_exported"]["satisfied"])
        matrix = self.payload["constructed_theoretical_objects"]["finite_acceptance_matrix"]
        self.assertEqual(len(matrix), 64)
        self.assertEqual(sum(1 for row in matrix if row["accepts_strict_p2938_provenance_theorem"]), 1)

    def test_certificate_lay_summary_nonpromotion_and_docs(self):
        cert = self.payload["provenance_certificate"]
        self.assertTrue(cert["K_plus_C_equals_V"])
        self.assertFalse(cert["strict_nadsoliton_functor_exported"])
        self.assertFalse(cert["nonconventional_equal_summand_law_exported"])
        self.assertFalse(cert["aggregate_localizer_exported"])
        self.assertFalse(cert["ratio_package_beta_unit_coupling_exported"])
        self.assertFalse(cert["p2951_p2938_strict_provenance_atom_discharged"])
        self.assertIn("real positive progress", self.payload["lay_summary"]["positive_progress"])
        self.assertIn("not yet a breakthrough", self.payload["lay_summary"]["not_a_breakthrough"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2958/S1908", MD.read_text(encoding="utf-8"))
        self.assertIn("P2958/S1908", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2958/S1908", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2958", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
