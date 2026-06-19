import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2938_s1888_unit_character_enriched_prime_coordinate_source_candidate.py"
OUT = ROOT / "generated" / "p2938_s1888_unit_character_enriched_prime_coordinate_source_candidate.json"
MD = ROOT / "generated" / "p2938_s1888_unit_character_enriched_prime_coordinate_source_candidate.md"


class P2938UnitCharacterEnrichedPrimeCoordinateSourceCandidateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2938_UNIT_CHARACTER_ENRICHED_PRIME_COORDINATE_SOURCE_CANDIDATE_CONDITIONAL_REJECTED")
        self.assertIsNotNone(self.payload["input_hashes"]["P2937"])

    def test_character_table_and_candidate_vector(self):
        cert = self.payload["candidate_certificate"]
        self.assertEqual(cert["unit_group"], [1, 5, 7, 11])
        self.assertEqual(cert["character_count"], 4)
        self.assertEqual(cert["nontrivial_character_count"], 3)
        self.assertEqual(cert["prime_coordinate_vector_order_2_3_5_7_11"], [1, 2, 2, 2, 2])
        self.assertTrue(cert["all_prime_coordinates_nonzero"])
        self.assertEqual(cert["product_pair_count_de_le_11"], 29)
        self.assertEqual(cert["product_additivity_defect_count"], 0)
        self.assertFalse(cert["accepted_strict_prime_log_source"])

    def test_rows_and_no_closure_flags(self):
        char_rows = self.payload["constructed_theoretical_objects"]["unit_character_rows"]
        self.assertEqual(len(char_rows), 4)
        self.assertTrue(all(row["multiplicative_on_U12"] for row in char_rows))
        prime_rows = self.payload["constructed_theoretical_objects"]["prime_coordinate_rows"]
        self.assertTrue(all(row["combined_coordinate"] != 0 for row in prime_rows))
        acceptance = self.payload["constructed_theoretical_objects"]["acceptance_rows"]
        self.assertFalse(next(row for row in acceptance if row["criterion"] == "strict_nadsoliton_formula_provenance")["satisfied"])
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_aut_breaking_prime_coordinate_source_law_exported", "strict_prime_log_value_source_exported", "strict_delta_eta_source_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2938/S1888", MD.read_text(encoding="utf-8"))
        self.assertIn("P2938/S1888", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2938/S1888", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2938", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
