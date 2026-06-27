import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3126_s2076_pi_point_pointed_support_source_audit.py"
OUT = ROOT / "generated" / "p3126_s2076_pi_point_pointed_support_source_audit.json"
MD = ROOT / "generated" / "p3126_s2076_pi_point_pointed_support_source_audit.md"


class P3126PiPointPointedSupportSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3126_PI_POINT_POINTED_SUPPORT_SOURCE_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3125"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3125_accepted_Lambda_origin_sources"], 0)
        self.assertEqual(cert["candidate_Pi_point_sources"], 19)
        self.assertEqual(cert["support_law_rows"], 323)
        self.assertEqual(cert["symmetry_witness_rows"], 228)
        self.assertEqual(cert["Lambda_Phi_Delta_Iota_Kappa_Tau_Xi_R_coupling_rows"], 209)
        self.assertEqual(cert["candidate_gate_rows"], 437)
        self.assertGreaterEqual(cert["promising_internal_pointed_supports"], 10)
        self.assertEqual(cert["accepted_Pi_point_sources"], 0)

    def test_candidate_matrix_and_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_Pi_point_sources"]
        self.assertTrue(any(row["candidate"] == "minimal_phase_information_support_point" and row["unique_support_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "minimal_phase_information_support_point" and not row["translation_safe"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "a_phi_balanced_support_point" and row["C_phi_A_phi_preserved"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "translation_orbit_support_class" and row["translation_safe"] and not row["unique_support_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "observed_light_event_support" and not row["strict_nadsoliton_data_only"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_chosen_support_point" and not row["selector_free"] for row in candidates))
        self.assertTrue(all(not row["accepted_Pi_point_source"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["accepted_coupling_chain"] for row in objs["Lambda_Phi_Delta_Iota_Kappa_Tau_Xi_R_coupling_rows"]))

    def test_recommendation_and_docs(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["internal_phase_information_pointed_supports_remain_promising_but_scoped"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("Internal phase/information support objects", decision["bounded_result"])
        self.assertIn("strict orbit-tie-breaking invariant object Omega_tie", decision["next_honest_step"])
        self.assertIn("without becoming a selector premise", decision["next_honest_step"])
        self.assertIn("P3126/S2076", MD.read_text(encoding="utf-8"))
        self.assertIn("P3126/S2076", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3126/S2076", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3126", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
