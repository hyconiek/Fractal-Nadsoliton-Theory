import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2955_s1905_p2601_identity_action_to_delta_numerator_source_map_obstruction.py"
OUT = ROOT / "generated" / "p2955_s1905_p2601_identity_action_to_delta_numerator_source_map_obstruction.json"
MD = ROOT / "generated" / "p2955_s1905_p2601_identity_action_to_delta_numerator_source_map_obstruction.md"


class P2955P2601IdentityActionToDeltaNumeratorSourceMapObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2955_P2601_IDENTITY_ACTION_TO_DELTA_NUMERATOR_SOURCE_MAP_OBSTRUCTION_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2601"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2954"])

    def test_source_map_certificate(self):
        cert = self.payload["source_map_certificate"]
        self.assertTrue(cert["p2601_identity_action_source_exported"])
        self.assertTrue(cert["p2954_target_signature_available"])
        self.assertEqual(cert["accepted_unital_multiplicative_row_count"], 3)
        self.assertFalse(cert["explicit_functor_exported"])
        self.assertFalse(cert["carrier_zero_alias_excluded_by_source_map"])
        self.assertFalse(cert["ratio_package_delta_numerator_source_exported"])
        self.assertFalse(cert["p2951_identity_deficit_atom_discharged"])
        self.assertEqual(cert["interface_matrix_row_count"], 4)
        self.assertEqual(cert["accepted_interface_matrix_row_count"], 0)

    def test_acceptance_rows_and_matrix(self):
        obj = self.payload["constructed_theoretical_objects"]
        summary = obj["p2601_export_summary"]
        self.assertEqual(summary["accepted_slopes"], ["1/2", "4/5", "3/2"])
        self.assertEqual(summary["accepted_intercepts"], ["0"])
        checks = {row["obligation"]: row for row in obj["source_map_acceptance_rows"]}
        self.assertTrue(checks["upstream_identity_action_source_exists"]["satisfied"])
        self.assertTrue(checks["target_identity_signature_available"]["satisfied"])
        self.assertFalse(checks["explicit_functor_from_p2601_source_to_p2954_delta_numerator_signature"]["satisfied"])
        self.assertFalse(checks["carrier_zero_alias_exclusion_under_source_map"]["satisfied"])
        self.assertFalse(checks["ratio_package_delta_numerator_source_exported"]["satisfied"])
        matrix = obj["finite_interface_matrix"]
        self.assertEqual(len(matrix), 4)
        self.assertFalse(any(row["accepts_strict_identity_deficit_source_map"] for row in matrix))

    def test_nonpromotion_and_docs_updated(self):
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["explicit_p2601_to_p2954_delta_numerator_functor_exported", "strict_identity_deficit_source_law_exported", "strict_delta_eta_source_law_exported", "strict_ratio_package_source_theorem_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])
        self.assertIn("P2955/S1905", MD.read_text(encoding="utf-8"))
        self.assertIn("P2955/S1905", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2955/S1905", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2955", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
