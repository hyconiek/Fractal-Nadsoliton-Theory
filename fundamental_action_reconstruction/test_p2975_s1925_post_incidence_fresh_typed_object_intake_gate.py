import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p2975_s1925_post_incidence_fresh_typed_object_intake_gate.py"
OUT = ROOT / "generated" / "p2975_s1925_post_incidence_fresh_typed_object_intake_gate.json"
MD = ROOT / "generated" / "p2975_s1925_post_incidence_fresh_typed_object_intake_gate.md"

class P2975PostIncidenceFreshTypedObjectIntakeGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P2975_POST_INCIDENCE_FRESH_TYPED_OBJECT_INTAKE_GATE_NO_NEW_LIVE_FRONTIER")
        self.assertIsNotNone(self.payload["input_hashes"]["P2974"])

    def test_intake_certificate(self):
        cert = self.payload["intake_certificate"]
        self.assertEqual(cert["candidate_count"], 7)
        self.assertEqual(cert["accepted_candidate_count"], 0)
        self.assertEqual(cert["lane_count"], 6)
        self.assertEqual(cert["unlocked_lane_count"], 0)
        self.assertEqual(cert["acceptance_matrix_rows"], 64)
        self.assertEqual(cert["accepted_matrix_rows"], 1)
        self.assertTrue(cert["no_new_live_frontier_certificate_exported"])

    def test_candidates_and_lanes(self):
        rows = {r["candidate"]: r for r in self.payload["constructed_theoretical_objects"]["candidate_intake_rows"]}
        self.assertFalse(rows["incidence_formal_derivative_table_replay"]["genuinely_new_outside_incidence_lane"])
        self.assertFalse(rows["ratio_package_scalar_9_5_or_unit_replay"]["not_replay_of_closed_ratio_or_gamma_or_selector_or_bridge_lane"])
        self.assertFalse(rows["fresh_strict_typed_object_placeholder"]["accepted_current_new_live_frontier"])
        lanes = {r["lane"]: r for r in self.payload["constructed_theoretical_objects"]["state_map_lane_rows"]}
        self.assertFalse(lanes["incidence theorem triad"]["live_frontier_unlocked"])
        self.assertFalse(lanes["L_total/ToE promotion lane"]["live_frontier_unlocked"])

    def test_docs_and_nonpromotion(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertIn("P2975/S1925", MD.read_text(encoding="utf-8"))
        self.assertIn("P2975/S1925", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2975/S1925", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P2975", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
