import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2950_s1900_exact_package_beta_eta_scale_coupling_obstruction.py"
OUT = ROOT / "generated" / "p2950_s1900_exact_package_beta_eta_scale_coupling_obstruction.json"
MD = ROOT / "generated" / "p2950_s1900_exact_package_beta_eta_scale_coupling_obstruction.md"


class P2950ExactPackageBetaEtaScaleCouplingObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2950_EXACT_PACKAGE_BETA_ETA_SCALE_COUPLING_OBSTRUCTION_NO_STRICT_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P2928"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2948"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2949"])

    def test_beta_eta_scale_coupling_certificate(self):
        cert = self.payload["beta_eta_scale_coupling_certificate"]
        self.assertEqual(cert["delta"]["as_string"], "4/5")
        self.assertEqual(cert["eta"]["as_string"], "9/5")
        self.assertTrue(cert["eta_equals_one_plus_delta"])
        self.assertEqual(cert["formal_multiplicative_carrier_defect_count"], 0)
        self.assertEqual(cert["sample_positive_beta_scale_count"], 4)
        self.assertTrue(cert["all_sample_beta_scales_leave_eta_package_unchanged"])
        self.assertFalse(cert["strict_positive_beta_scale_source_exported"])
        self.assertFalse(cert["p2948_beta_eta_coupling_premise_discharged"])

    def test_rows_acceptance_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertTrue(all(row["satisfied_finitely"] for row in obj["exact_package_rows"]))
        beta_rows = obj["positive_beta_scale_orbit_rows"]
        self.assertEqual({row["beta_scale"]["as_string"] for row in beta_rows}, {"1/2", "1/1", "2/1", "5/1"})
        self.assertTrue(all(row["same_eta_ratio_package"] for row in beta_rows))
        self.assertTrue(all(not row["strict_source_selects_this_beta"] for row in beta_rows))
        obligations = obj["theorem_obligation_rows"]
        self.assertTrue(obligations[0]["satisfied"])
        self.assertTrue(obligations[1]["satisfied"])
        self.assertFalse(any(row["satisfied"] for row in obligations[2:]))
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_beta_eta_coupling_theorem_exported", "strict_positive_beta_scale_source_exported", "strict_damping_beta_eta_source_packet_exported", "strict_delta_eta_source_law_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2950/S1900", MD.read_text(encoding="utf-8"))
        self.assertIn("P2950/S1900", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2950/S1900", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2950", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
