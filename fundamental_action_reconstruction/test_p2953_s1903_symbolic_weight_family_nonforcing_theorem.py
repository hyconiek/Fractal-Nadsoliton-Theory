import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2953_s1903_symbolic_weight_family_nonforcing_theorem.py"
OUT = ROOT / "generated" / "p2953_s1903_symbolic_weight_family_nonforcing_theorem.json"
MD = ROOT / "generated" / "p2953_s1903_symbolic_weight_family_nonforcing_theorem.md"


class P2953SymbolicWeightFamilyNonforcingTheoremTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2953_SYMBOLIC_WEIGHT_FAMILY_NONFORCING_THEOREM_NO_STRICT_SOURCE")
        self.assertIsNotNone(self.payload["input_hashes"]["P2938"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2952"])

    def test_symbolic_weight_certificate(self):
        cert = self.payload["symbolic_weight_certificate"]
        self.assertEqual(cert["closed_form"], "V(a,b)=[a,2a,2b,2b,2b]")
        self.assertTrue(cert["unbounded_positive_rational_family_constructed"])
        self.assertEqual(cert["sample_non_target_positive_witness_count"], 3)
        self.assertTrue(cert["target_vector_equivalent_to_a_equals_1_b_equals_1"])
        self.assertTrue(cert["p2952_bounded_scan_caveat_removed"])
        self.assertFalse(cert["strict_equal_weight_source_theorem_exported"])
        self.assertFalse(cert["p2951_torsion_character_provenance_atom_discharged"])

    def test_symbolic_rows_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        vectors = obj["ingredient_vectors"]
        self.assertEqual(vectors["kernel_excess_vector_K_order_2_3_5_7_11"], [1, 2, 0, 0, 0])
        self.assertEqual(vectors["character_negativity_vector_C_order_2_3_5_7_11"], [0, 0, 2, 2, 2])
        theorem_rows = obj["symbolic_theorem_rows"]
        self.assertTrue(all(row["verified"] for row in theorem_rows[:4]))
        self.assertFalse(theorem_rows[4]["verified"])
        witnesses = obj["positive_rational_witness_rows"]
        self.assertTrue(all(row["positive_coordinates"] for row in witnesses))
        self.assertFalse(any(row["matches_p2938_target"] for row in witnesses))
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_equal_weight_source_theorem_exported", "strict_torsion_character_source_theorem_exported", "strict_ratio_package_source_theorem_exported", "strict_delta_eta_source_law_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2953/S1903", MD.read_text(encoding="utf-8"))
        self.assertIn("P2953/S1903", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2953/S1903", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2953", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
