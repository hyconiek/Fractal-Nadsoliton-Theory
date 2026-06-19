import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2929_s1879_post_damping_state_map_no_new_live_frontier_certificate.py"
OUT = ROOT / "generated" / "p2929_s1879_post_damping_state_map_no_new_live_frontier_certificate.json"
MD = ROOT / "generated" / "p2929_s1879_post_damping_state_map_no_new_live_frontier_certificate.md"


class P2929PostDampingStateMapNoNewLiveFrontierCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P2929_POST_DAMPING_STATE_MAP_NO_NEW_LIVE_FRONTIER_CERTIFICATE")
        self.assertIsNotNone(self.payload["input_hashes"]["P2920"])
        self.assertIsNotNone(self.payload["input_hashes"]["P2928"])

    def test_state_map_counts(self):
        cert = self.payload["state_map_certificate"]
        self.assertEqual(cert["lane_count"], 9)
        self.assertEqual(cert["live_frontier_unlocked_lane_count"], 0)
        self.assertEqual(cert["candidate_intake_count"], 5)
        self.assertEqual(cert["accepted_new_live_frontier_candidate_count"], 0)
        self.assertTrue(cert["all_lanes_preserve_nonpromotion"])
        self.assertTrue(cert["no_new_live_frontier_certificate_exported"])

    def test_intake_gate_and_false_closure_flags(self):
        gate = self.payload["constructed_theoretical_objects"]["fresh_typed_object_intake_gate"]
        self.assertEqual(gate["name"], "Fresh_Strict_Typed_Object_Intake_Gate")
        self.assertEqual(len(gate["obligations"]), 5)
        flags = self.payload["decision"]["negative_export_flags"]
        for key in ["new_live_frontier_unlocked", "strict_damping_beta_eta_source_packet_exported", "strict_selector_closure_exported", "bridge_closure_exported", "role_transfer_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "toe_closure_exported"]:
            self.assertFalse(flags[key])

    def test_docs_updated(self):
        self.assertIn("P2929/S1879", MD.read_text(encoding="utf-8"))
        self.assertIn("P2929/S1879", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2929/S1879", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2929", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
