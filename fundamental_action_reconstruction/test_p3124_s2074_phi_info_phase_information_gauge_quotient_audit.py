import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3124_s2074_phi_info_phase_information_gauge_quotient_audit.py"
OUT = ROOT / "generated" / "p3124_s2074_phi_info_phase_information_gauge_quotient_audit.json"
MD = ROOT / "generated" / "p3124_s2074_phi_info_phase_information_gauge_quotient_audit.md"


class P3124PhiInfoPhaseInformationGaugeQuotientAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3124_PHI_INFO_PHASE_INFORMATION_GAUGE_QUOTIENT_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3123"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3123_accepted_Delta_asym_sources"], 0)
        self.assertEqual(cert["candidate_Phi_Info_sources"], 17)
        self.assertEqual(cert["quotient_law_rows"], 221)
        self.assertEqual(cert["gauge_orbit_rows"], 136)
        self.assertEqual(cert["Delta_Iota_Kappa_Tau_Xi_R_coupling_rows"], 153)
        self.assertEqual(cert["candidate_gate_rows"], 323)
        self.assertGreaterEqual(cert["internal_promising_phase_info_candidates"], 5)
        self.assertEqual(cert["accepted_Phi_Info_sources"], 0)

    def test_candidate_matrix_and_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_Phi_Info_sources"]
        self.assertTrue(any(row["candidate"] == "phase_information_orbit_quotient" and row["phase_origin_gauge_quotiented"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "phase_information_orbit_quotient" and not row["Delta_asym_source_unlocked"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "a_phi_normalized_info_gradient" and row["C_phi_A_phi_preserved"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "observed_light_phase_anchor" and not row["strict_nadsoliton_data_only"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_fixed_phase_origin" and not row["selector_bridge_ltotal_toe_free"] for row in candidates))
        self.assertTrue(all(not row["accepted_Phi_Info_source"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["accepted_coupling_chain"] for row in objs["Delta_Iota_Kappa_Tau_Xi_R_coupling_rows"]))

    def test_recommendation_and_docs(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["phase_information_lane_confirmed_as_best_scoped_shape"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("phase-information remains the strongest scoped internal lane", decision["bounded_result"])
        self.assertIn("strict phase-origin/source-localizer object Lambda_origin", decision["next_honest_step"])
        self.assertIn("nonzero phase-information quotient representative", decision["next_honest_step"])
        self.assertIn("P3124/S2074", MD.read_text(encoding="utf-8"))
        self.assertIn("P3124/S2074", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3124/S2074", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3124", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
