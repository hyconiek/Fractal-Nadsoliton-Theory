import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2952_s1902_torsion_character_aggregate_weight_source_obstruction.py"
OUT = ROOT / "generated" / "p2952_s1902_torsion_character_aggregate_weight_source_obstruction.json"
MD = ROOT / "generated" / "p2952_s1902_torsion_character_aggregate_weight_source_obstruction.md"


class P2952TorsionCharacterAggregateWeightSourceObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2952_TORSION_CHARACTER_AGGREGATE_WEIGHT_SOURCE_OBSTRUCTION_NO_STRICT_PROVENANCE")
        self.assertIsNotNone(self.payload["input_hashes"]["P2938"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2951"])

    def test_weight_source_certificate(self):
        cert = self.payload["weight_source_certificate"]
        self.assertEqual(cert["weight_domain"], [1, 2, 3, 4])
        self.assertEqual(cert["family_row_count"], 16)
        self.assertTrue(cert["all_family_rows_product_additive_by_construction"])
        self.assertEqual(cert["target_weight_pair_count_in_domain"], 1)
        self.assertEqual(cert["target_weight_pairs"], [[1, 1]])
        self.assertTrue(cert["target_vector_forces_equal_weights_in_a_b_family"])
        self.assertFalse(cert["strict_equal_weight_source_theorem_exported"])
        self.assertFalse(cert["p2951_torsion_character_provenance_atom_discharged"])

    def test_rows_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual([row["p2938_equal_weight_coordinate"] for row in obj["ingredient_rows"]], [1, 2, 2, 2, 2])
        family = obj["bounded_positive_weight_family_rows"]
        self.assertEqual(len([row for row in family if row["matches_p2938_target_vector"]]), 1)
        self.assertTrue(all(row["all_coordinates_positive"] for row in family))
        obligations = obj["theorem_obligation_rows"]
        self.assertTrue(obligations[0]["satisfied"])
        self.assertTrue(obligations[1]["satisfied"])
        self.assertFalse(any(row["satisfied"] for row in obligations[2:]))
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_equal_weight_source_theorem_exported", "strict_torsion_character_source_theorem_exported", "strict_ratio_package_source_theorem_exported", "strict_delta_eta_source_law_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2952/S1902", MD.read_text(encoding="utf-8"))
        self.assertIn("P2952/S1902", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2952/S1902", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2952", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
