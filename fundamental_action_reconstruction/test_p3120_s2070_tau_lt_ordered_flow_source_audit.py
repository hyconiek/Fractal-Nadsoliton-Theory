import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3120_s2070_tau_lt_ordered_flow_source_audit.py"
OUT = ROOT / "generated" / "p3120_s2070_tau_lt_ordered_flow_source_audit.json"
MD = ROOT / "generated" / "p3120_s2070_tau_lt_ordered_flow_source_audit.md"


class P3120TauLTOrderedFlowSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3120_TAU_LT_ORDERED_FLOW_SOURCE_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3119"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3119_accepted_Xi_LT_axis_sources"], 0)
        self.assertEqual(cert["candidate_Tau_LT_ordered_flows"], 13)
        self.assertEqual(cert["flow_law_rows"], 117)
        self.assertEqual(cert["flow_scale_covariance_rows"], 65)
        self.assertEqual(cert["Xi_LT_R_dim_coupling_rows"], 91)
        self.assertEqual(cert["candidate_gate_rows"], 195)
        self.assertEqual(cert["accepted_Tau_LT_ordered_flows"], 0)

    def test_candidate_matrix_and_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_Tau_LT_ordered_flows"]
        self.assertTrue(any(row["candidate"] == "entropy_monotone_flow_bifunctor" and row["U_time_as_internal_order_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "phase_winding_order_extension" and not row["cycle_or_recurrence_cut_nonconventional"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "lightcone_flow_bifunctor" and not row["strict_nadsoliton_data_only"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_oriented_flow_bifunctor" and not row["selector_bridge_ltotal_toe_free"] for row in candidates))
        self.assertTrue(all(not row["accepted_Tau_LT_ordered_flow"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["accepted_coupling_chain"] for row in objs["Xi_LT_R_dim_coupling_rows"]))

    def test_negative_flags_docs_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["candidate_Tau_LT_ordered_flows_constructed"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("strict recurrence-cut source object Kappa_cycle", decision["next_honest_step"])
        self.assertIn("acyclicity/irreversibility law", decision["next_honest_step"])
        self.assertIn("P3120/S2070", MD.read_text(encoding="utf-8"))
        self.assertIn("P3120/S2070", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3120/S2070", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3120", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
