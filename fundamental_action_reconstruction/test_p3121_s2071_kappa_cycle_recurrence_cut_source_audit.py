import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3121_s2071_kappa_cycle_recurrence_cut_source_audit.py"
OUT = ROOT / "generated" / "p3121_s2071_kappa_cycle_recurrence_cut_source_audit.json"
MD = ROOT / "generated" / "p3121_s2071_kappa_cycle_recurrence_cut_source_audit.md"


class P3121KappaCycleRecurrenceCutSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3121_KAPPA_CYCLE_RECURRENCE_CUT_SOURCE_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3120"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3120_accepted_Tau_LT_ordered_flows"], 0)
        self.assertEqual(cert["candidate_Kappa_cycle_sources"], 14)
        self.assertEqual(cert["cycle_law_rows"], 140)
        self.assertEqual(cert["finite_graph_witness_rows"], 84)
        self.assertEqual(cert["Tau_LT_Xi_LT_R_dim_coupling_rows"], 112)
        self.assertEqual(cert["candidate_gate_rows"], 224)
        self.assertEqual(cert["accepted_Kappa_cycle_sources"], 0)

    def test_candidate_matrix_and_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_Kappa_cycle_sources"]
        self.assertTrue(any(row["candidate"] == "entropy_strict_increase_cut" and row["cycle_cut_law_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "z12_lift_acyclic_cover_cut" and not row["gauge_invariant_not_label_choice"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "thermodynamic_environment_cut" and not row["strict_nadsoliton_data_only"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_oriented_cut" and not row["selector_bridge_ltotal_toe_free"] for row in candidates))
        self.assertTrue(all(not row["accepted_Kappa_cycle_source"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["accepted_coupling_chain"] for row in objs["Tau_LT_Xi_LT_R_dim_coupling_rows"]))

    def test_negative_flags_docs_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["candidate_Kappa_cycle_sources_constructed"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("strict irreversible-defect source object Iota_irrev", decision["next_honest_step"])
        self.assertIn("nonzero signed defect functional", decision["next_honest_step"])
        self.assertIn("P3121/S2071", MD.read_text(encoding="utf-8"))
        self.assertIn("P3121/S2071", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3121/S2071", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3121", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
