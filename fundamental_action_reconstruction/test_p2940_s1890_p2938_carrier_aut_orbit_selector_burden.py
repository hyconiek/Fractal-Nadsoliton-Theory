import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2940_s1890_p2938_carrier_aut_orbit_selector_burden.py"
OUT = ROOT / "generated" / "p2940_s1890_p2938_carrier_aut_orbit_selector_burden.json"
MD = ROOT / "generated" / "p2940_s1890_p2938_carrier_aut_orbit_selector_burden.md"


class P2940P2938CarrierAutOrbitSelectorBurdenTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2940_P2938_CARRIER_AUT_ORBIT_SELECTOR_BURDEN_NO_STRICT_PROVENANCE")
        self.assertIsNotNone(self.payload["input_hashes"]["P2939"])

    def test_orbit_certificate(self):
        cert = self.payload["orbit_certificate"]
        self.assertEqual(cert["orbit_count"], 5)
        self.assertEqual(cert["selector_burden_orbit_count"], 4)
        self.assertEqual(cert["orbit_constant_count"], 1)
        self.assertFalse(cert["all_orbits_constant"])
        self.assertFalse(cert["strict_selector_source_exported"])
        self.assertFalse(cert["strict_provenance_theorem_exported"])

    def test_orbit_rows_and_flags(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["prime_coordinate_vector_order_2_3_5_7_11"], [1, 2, 2, 2, 2])
        burden_rows = obj["selector_burden_rows"]
        self.assertEqual(len(burden_rows), 4)
        self.assertTrue(any(row["orbit"] == [1, 5, 7, 11] for row in burden_rows))
        criteria = obj["acceptance_rows"]
        self.assertTrue(next(row for row in criteria if row["criterion"] == "selector_burden_quantified")["satisfied"])
        self.assertFalse(next(row for row in criteria if row["criterion"] == "strict_selector_or_symmetry_breaking_source_exported")["satisfied"])
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_selector_source_exported", "strict_aut_breaking_prime_coordinate_source_law_exported", "strict_prime_log_value_source_exported", "strict_delta_eta_source_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2940/S1890", MD.read_text(encoding="utf-8"))
        self.assertIn("P2940/S1890", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2940/S1890", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2940", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
