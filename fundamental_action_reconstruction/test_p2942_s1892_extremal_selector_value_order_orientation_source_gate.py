import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2942_s1892_extremal_selector_value_order_orientation_source_gate.py"
OUT = ROOT / "generated" / "p2942_s1892_extremal_selector_value_order_orientation_source_gate.json"
MD = ROOT / "generated" / "p2942_s1892_extremal_selector_value_order_orientation_source_gate.md"


class P2942ExtremalSelectorValueOrderOrientationSourceGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2942_EXTREMAL_SELECTOR_VALUE_ORDER_ORIENTATION_SOURCE_GATE_NO_STRICT_POLARITY_SOURCE")
        self.assertIsNotNone(self.payload["input_hashes"]["P2941"])

    def test_orientation_certificate(self):
        cert = self.payload["orientation_certificate"]
        self.assertEqual(cert["scan_row_count"], 35)
        self.assertEqual(cert["positive_orientation_row_count"], 15)
        self.assertEqual(cert["negative_orientation_row_count"], 15)
        self.assertEqual(cert["zero_orientation_row_count"], 5)
        self.assertEqual(cert["positive_rows_unique_on_all_orbits"], 15)
        self.assertEqual(cert["negative_rows_unique_on_all_orbits"], 0)
        self.assertEqual(cert["zero_rows_unique_on_all_orbits"], 0)
        self.assertFalse(cert["strict_positive_orientation_source_theorem_exported"])
        self.assertFalse(cert["strict_provenance_theorem_exported"])
        self.assertFalse(cert["accepted_strict_source"])

    def test_acceptance_split_and_nonpromotion(self):
        obj = self.payload["constructed_theoretical_objects"]
        criteria = obj["acceptance_rows"]
        self.assertTrue(next(row for row in criteria if row["criterion"] == "affine_value_order_orientation_space_computed")["satisfied"])
        self.assertTrue(next(row for row in criteria if row["criterion"] == "positive_orientation_conditionally_selects_unique_min_skeleton")["satisfied"])
        self.assertTrue(next(row for row in criteria if row["criterion"] == "negative_orientation_rejected_by_unit_orbit_tie")["satisfied"])
        self.assertTrue(next(row for row in criteria if row["criterion"] == "zero_orientation_rejected_by_collapse")["satisfied"])
        self.assertFalse(next(row for row in criteria if row["criterion"] == "strict_positive_orientation_source_theorem_exported")["satisfied"])
        self.assertFalse(next(row for row in criteria if row["criterion"] == "delta_eta_and_beta_eta_coupling_exported")["satisfied"])
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["strict_positive_orientation_source_theorem_exported", "strict_selector_source_exported", "strict_aut_breaking_prime_coordinate_source_law_exported", "strict_prime_log_value_source_exported", "strict_delta_eta_source_exported", "strict_beta_eta_coupling_theorem_exported", "strict_damping_beta_eta_source_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_candidate_rows_and_docs_updated(self):
        candidates = self.payload["constructed_theoretical_objects"]["source_candidate_rows"]
        self.assertEqual(len(candidates), 5)
        self.assertTrue(any(row["candidate"] == "positive_affine_value_order_q(v)=a*v+b_with_a>0" for row in candidates))
        self.assertIn("P2942/S1892", MD.read_text(encoding="utf-8"))
        self.assertIn("P2942/S1892", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2942/S1892", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2942", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
