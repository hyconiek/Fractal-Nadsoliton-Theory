import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2939_s1889_p2938_carrier_aut_breaking_provenance_boundary.py"
OUT = ROOT / "generated" / "p2939_s1889_p2938_carrier_aut_breaking_provenance_boundary.json"
MD = ROOT / "generated" / "p2939_s1889_p2938_carrier_aut_breaking_provenance_boundary.md"


class P2939P2938CarrierAutBreakingProvenanceBoundaryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2939_P2938_CARRIER_AUT_BREAKING_PROVENANCE_BOUNDARY_NO_STRICT_SOURCE")
        self.assertIsNotNone(self.payload["input_hashes"]["P2938"])

    def test_boundary_certificate(self):
        cert = self.payload["boundary_certificate"]
        self.assertEqual(cert["product_pair_count_de_le_11"], 29)
        self.assertEqual(cert["product_additivity_defect_count"], 0)
        self.assertEqual(cert["aut_action_row_count"], 44)
        self.assertGreater(cert["nontrivial_aut_defect_count"], 0)
        self.assertIsNotNone(cert["first_aut_breaking_witness"])
        self.assertEqual(cert["algebraic_readiness_criteria_satisfied_count"], 4)
        self.assertEqual(cert["strict_theorem_criteria_satisfied_count"], 0)
        self.assertFalse(cert["accepted_strict_source"])

    def test_aut_breaking_and_no_closure_flags(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["prime_coordinate_vector_order_2_3_5_7_11"], [1, 2, 2, 2, 2])
        self.assertTrue(any(not row["preserves_value"] for row in obj["aut_action_defect_rows"] if row["unit"] != 1))
        criteria = obj["provenance_acceptance_rows"]
        self.assertTrue(next(row for row in criteria if row["criterion"] == "AutZ12_breaking_witness_exists")["satisfied"])
        self.assertFalse(next(row for row in criteria if row["criterion"] == "strict_nadsoliton_provenance_theorem")["satisfied"])
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_aut_breaking_prime_coordinate_source_law_exported", "strict_prime_log_value_source_exported", "strict_delta_eta_source_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2939/S1889", MD.read_text(encoding="utf-8"))
        self.assertIn("P2939/S1889", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2939/S1889", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2939", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
