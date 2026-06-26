import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3119_s2069_xi_lt_axis_source_object_audit.py"
OUT = ROOT / "generated" / "p3119_s2069_xi_lt_axis_source_object_audit.json"
MD = ROOT / "generated" / "p3119_s2069_xi_lt_axis_source_object_audit.md"


class P3119XiLTAxisSourceObjectAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3119_XI_LT_AXIS_SOURCE_OBJECT_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3118"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3118_accepted_R_dim_relations"], 0)
        self.assertEqual(cert["candidate_Xi_LT_axis_sources"], 12)
        self.assertEqual(cert["axis_source_rows"], 96)
        self.assertEqual(cert["axis_scale_covariance_rows"], 60)
        self.assertEqual(cert["R_dim_coupling_rows"], 72)
        self.assertEqual(cert["candidate_gate_rows"], 168)
        self.assertEqual(cert["accepted_Xi_LT_axis_sources"], 0)

    def test_candidate_matrix_and_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_Xi_LT_axis_sources"]
        self.assertTrue(any(row["candidate"] == "phase_gradient_axis_split" and row["U_length_source_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "symplectic_conjugate_axis_pair" and not row["not_unit_convention"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "planck_light_axis_pair" and not row["strict_nadsoliton_data_only"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_oriented_axis_pair" and not row["selector_bridge_ltotal_toe_free"] for row in candidates))
        self.assertTrue(all(not row["accepted_Xi_LT_axis_source"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["accepted_coupling_chain"] for row in objs["R_dim_coupling_rows"]))

    def test_negative_flags_docs_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["candidate_Xi_LT_axis_sources_constructed"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("strict ordered-flow source object Tau_LT", decision["next_honest_step"])
        self.assertIn("temporal-order/length-extension bifunctor", decision["next_honest_step"])
        self.assertIn("P3119/S2069", MD.read_text(encoding="utf-8"))
        self.assertIn("P3119/S2069", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3119/S2069", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3119", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
