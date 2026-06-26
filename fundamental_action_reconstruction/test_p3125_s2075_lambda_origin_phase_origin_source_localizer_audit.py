import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3125_s2075_lambda_origin_phase_origin_source_localizer_audit.py"
OUT = ROOT / "generated" / "p3125_s2075_lambda_origin_phase_origin_source_localizer_audit.json"
MD = ROOT / "generated" / "p3125_s2075_lambda_origin_phase_origin_source_localizer_audit.md"


class P3125LambdaOriginPhaseOriginSourceLocalizerAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3125_LAMBDA_ORIGIN_PHASE_ORIGIN_SOURCE_LOCALIZER_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3124"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3124_accepted_Phi_Info_sources"], 0)
        self.assertEqual(cert["candidate_Lambda_origin_sources"], 18)
        self.assertEqual(cert["localizer_law_rows"], 270)
        self.assertEqual(cert["symmetry_witness_rows"], 162)
        self.assertEqual(cert["Phi_Delta_Iota_Kappa_Tau_Xi_R_coupling_rows"], 180)
        self.assertEqual(cert["candidate_gate_rows"], 378)
        self.assertGreaterEqual(cert["promising_internal_localizers"], 8)
        self.assertEqual(cert["accepted_Lambda_origin_sources"], 0)

    def test_candidate_matrix_and_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_Lambda_origin_sources"]
        self.assertTrue(any(row["candidate"] == "phase_information_cross_extremum" and row["nonzero_representative_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "phase_information_cross_extremum" and not row["translation_gauge_safe"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "a_phi_cell_boundary_anchor" and row["C_phi_A_phi_preserved"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "observed_light_phase_event_anchor" and not row["strict_nadsoliton_data_only"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_chosen_origin_anchor" and not row["selector_bridge_ltotal_toe_free"] for row in candidates))
        self.assertTrue(all(not row["accepted_Lambda_origin_source"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["accepted_coupling_chain"] for row in objs["Phi_Delta_Iota_Kappa_Tau_Xi_R_coupling_rows"]))

    def test_recommendation_and_docs(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["internal_phase_information_localizers_remain_promising_but_scoped"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("Internal phase/information localizers", decision["bounded_result"])
        self.assertIn("strict pointed-support source object Pi_point", decision["next_honest_step"])
        self.assertIn("translation/inversion-safe pointed support theorem", decision["next_honest_step"])
        self.assertIn("P3125/S2075", MD.read_text(encoding="utf-8"))
        self.assertIn("P3125/S2075", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3125/S2075", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3125", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
