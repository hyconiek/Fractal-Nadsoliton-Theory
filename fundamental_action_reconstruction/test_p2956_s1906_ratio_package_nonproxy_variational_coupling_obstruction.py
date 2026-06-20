import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2956_s1906_ratio_package_nonproxy_variational_coupling_obstruction.py"
OUT = ROOT / "generated" / "p2956_s1906_ratio_package_nonproxy_variational_coupling_obstruction.json"
MD = ROOT / "generated" / "p2956_s1906_ratio_package_nonproxy_variational_coupling_obstruction.md"


class P2956RatioPackageNonproxyVariationalCouplingObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2956_RATIO_PACKAGE_NONPROXY_VARIATIONAL_COUPLING_OBSTRUCTION_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2948"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2951"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2955"])

    def test_variational_coupling_certificate(self):
        cert = self.payload["variational_coupling_certificate"]
        self.assertTrue(cert["finite_ratio_package_available"])
        self.assertTrue(cert["p2951_nonproxy_variational_atom_named"])
        self.assertFalse(cert["identity_deficit_source_map_closed"])
        self.assertTrue(cert["euler_rows_force_eta_only_by_9_to_5_insertion"])
        self.assertFalse(cert["independent_nonproxy_field_variable_exported"])
        self.assertFalse(cert["unit_bearing_action_density_exported"])
        self.assertFalse(cert["euler_forces_eta_without_insertion"])
        self.assertFalse(cert["p2951_nonproxy_variational_atom_discharged"])
        self.assertEqual(cert["acceptance_matrix_row_count"], 16)
        self.assertEqual(cert["acceptance_matrix_accepted_row_count"], 1)

    def test_package_euler_rows_and_matrix(self):
        obj = self.payload["constructed_theoretical_objects"]
        package = {row["object"]: row for row in obj["exact_package_rows"]}
        self.assertEqual(package["delta"]["value"]["as_string"], "4/5")
        self.assertEqual(package["eta"]["value"]["as_string"], "9/5")
        self.assertTrue(package["eta_equals_one_plus_delta"]["value"])
        euler_rows = obj["euler_source_ratio_rows"]
        self.assertEqual(len(euler_rows), 4)
        self.assertTrue(all(row["forces_eta_9_5"] for row in euler_rows))
        self.assertTrue(all(row["imports_9_to_5_source_ratio"] for row in euler_rows))
        self.assertFalse(any(row["accepted_as_strict_nonproxy_coupling"] for row in euler_rows))
        matrix = obj["finite_acceptance_matrix"]
        self.assertEqual(sum(1 for row in matrix if row["accepts_nonproxy_variational_damping_coupling"]), 1)
        self.assertEqual([row for row in matrix if row["accepts_nonproxy_variational_damping_coupling"]][0]["mask"], 15)

    def test_nonpromotion_and_docs_updated(self):
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_nonproxy_variational_damping_coupling_exported", "strict_ratio_package_source_theorem_exported", "strict_delta_eta_source_law_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])
        self.assertIn("P2956/S1906", MD.read_text(encoding="utf-8"))
        self.assertIn("P2956/S1906", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2956/S1906", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2956", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
