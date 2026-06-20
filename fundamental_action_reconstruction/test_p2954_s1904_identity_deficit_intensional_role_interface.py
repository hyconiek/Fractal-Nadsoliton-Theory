import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2954_s1904_identity_deficit_intensional_role_interface.py"
OUT = ROOT / "generated" / "p2954_s1904_identity_deficit_intensional_role_interface.json"
MD = ROOT / "generated" / "p2954_s1904_identity_deficit_intensional_role_interface.md"


class P2954IdentityDeficitIntensionalRoleInterfaceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2954_IDENTITY_DEFICIT_INTENSIONAL_ROLE_INTERFACE_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2938"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2944"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2949"])

    def test_intensional_role_certificate(self):
        cert = self.payload["intensional_role_certificate"]
        self.assertEqual(cert["node_count"], 11)
        self.assertEqual(cert["identity_extension"], [1])
        self.assertEqual(cert["carrier_zero_extension"], [1])
        self.assertTrue(cert["extensions_coincide"])
        self.assertTrue(cert["typed_role_signatures_separated"])
        self.assertTrue(cert["count_alias_replay_avoided"])
        self.assertFalse(cert["strict_identity_deficit_source_law_exported"])
        self.assertFalse(cert["p2951_identity_deficit_atom_discharged"])

    def test_interface_rows_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        checks = {row["interface_check"]: row for row in obj["interface_obligation_rows"]}
        self.assertTrue(checks["extension_coincidence_reproduced"]["satisfied"])
        self.assertTrue(checks["role_signatures_are_typed_distinct"]["satisfied"])
        self.assertTrue(checks["not_a_new_count_alias"]["satisfied"])
        self.assertFalse(checks["strict_identity_deficit_source_law_exported"]["satisfied"])
        first = obj["role_signature_rows"][0]
        self.assertEqual(first["node"], 1)
        self.assertEqual(first["identity_role_signature"]["defining_data_type"], "partial_monoid_action_law")
        self.assertEqual(first["carrier_zero_role_signature"]["defining_data_type"], "scalar_carrier_level_set")
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_identity_deficit_source_law_exported", "strict_delta_eta_source_law_exported", "strict_ratio_package_source_theorem_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2954/S1904", MD.read_text(encoding="utf-8"))
        self.assertIn("P2954/S1904", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2954/S1904", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2954", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
