import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2941_s1891_p2938_carrier_extremal_selector_polarity_obstruction.py"
OUT = ROOT / "generated" / "p2941_s1891_p2938_carrier_extremal_selector_polarity_obstruction.json"
MD = ROOT / "generated" / "p2941_s1891_p2938_carrier_extremal_selector_polarity_obstruction.md"


class P2941P2938CarrierExtremalSelectorPolarityObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2941_P2938_CARRIER_EXTREMAL_SELECTOR_POLARITY_OBSTRUCTION_NO_STRICT_PROVENANCE")
        self.assertIsNotNone(self.payload["input_hashes"]["P2940"])

    def test_selector_certificate(self):
        cert = self.payload["selector_certificate"]
        self.assertEqual(cert["orbit_selector_rows"], 10)
        self.assertTrue(cert["min_polarity_unique_on_all_orbits"])
        self.assertEqual(cert["max_polarity_tie_count"], 1)
        self.assertEqual(cert["u12_motion_rows"], 36)
        self.assertGreater(cert["u12_motion_defect_count"], 0)
        self.assertFalse(cert["strict_polarity_source_theorem_exported"])
        self.assertFalse(cert["strict_provenance_theorem_exported"])
        self.assertFalse(cert["accepted_strict_source"])

    def test_acceptance_and_nonpromotion_flags(self):
        obj = self.payload["constructed_theoretical_objects"]
        criteria = obj["acceptance_rows"]
        self.assertTrue(next(row for row in criteria if row["criterion"] == "carrier_extremal_selector_skeleton_constructed")["satisfied"])
        self.assertTrue(next(row for row in criteria if row["criterion"] == "min_polarity_unique_on_all_orbits")["satisfied"])
        self.assertFalse(next(row for row in criteria if row["criterion"] == "max_polarity_unique_on_all_orbits")["satisfied"])
        self.assertFalse(next(row for row in criteria if row["criterion"] == "global_section_U12_compatible")["satisfied"])
        self.assertFalse(next(row for row in criteria if row["criterion"] == "strict_polarity_source_theorem_exported")["satisfied"])
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_selector_source_exported", "strict_aut_breaking_prime_coordinate_source_law_exported", "strict_prime_log_value_source_exported", "strict_delta_eta_source_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2941/S1891", MD.read_text(encoding="utf-8"))
        self.assertIn("P2941/S1891", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2941/S1891", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2941", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
