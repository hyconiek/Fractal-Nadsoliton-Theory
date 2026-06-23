import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3043_s1993_post_selector_candidate_new_source_intake_gate.py"
OUT = ROOT / "generated" / "p3043_s1993_post_selector_candidate_new_source_intake_gate.json"
MD = ROOT / "generated" / "p3043_s1993_post_selector_candidate_new_source_intake_gate.md"

class P3043PostSelectorCandidateNewSourceIntakeGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_inputs(self):
        self.assertEqual(self.payload["status"], "P3043_POST_SELECTOR_CANDIDATE_NEW_SOURCE_INTAKE_NO_NEW_LIVE_FRONTIER")
        self.assertIsNotNone(self.payload["input_hashes"]["P3037"])
        self.assertIsNotNone(self.payload["input_hashes"]["P3042"])

    def test_intake_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["state_lanes"], 3)
        self.assertEqual(cert["exhausted_receiver_classes"], 5)
        self.assertEqual(cert["candidate_rows"], 4)
        self.assertEqual(cert["predicate_count"], 7)
        self.assertEqual(cert["accepted_new_source_law_rows"], 0)
        self.assertEqual(cert["replay_gated_rows"], 3)
        self.assertEqual(cert["unsupplied_new_source_slots"], 1)
        self.assertFalse(cert["new_live_frontier_unlocked"])

    def test_obligations_and_exports(self):
        obj = self.payload["constructed_theoretical_objects"]
        self.assertEqual(obj["object"], "PostSelectorCandidate_NewSourceIntakeGate")
        obligations = {row["obligation"]: row["satisfied"] for row in obj["proof_obligations"]}
        self.assertTrue(obligations["broad_state_map_consulted"])
        self.assertTrue(obligations["exhausted_receiver_classes_declared"])
        self.assertTrue(obligations["new_source_predicate_matrix_constructed"])
        self.assertFalse(obligations["concrete_new_strict_source_law_supplied"])
        self.assertFalse(obligations["selector_lane_reopened"])
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3043/S1993", MD.read_text(encoding="utf-8"))
        self.assertIn("P3043/S1993", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3043/S1993", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3043", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
