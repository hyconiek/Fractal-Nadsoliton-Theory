import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3118_s2068_r_dim_action_length_time_relation_audit.py"
OUT = ROOT / "generated" / "p3118_s2068_r_dim_action_length_time_relation_audit.json"
MD = ROOT / "generated" / "p3118_s2068_r_dim_action_length_time_relation_audit.md"


class P3118RDimActionLengthTimeRelationAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3118_R_DIM_ACTION_LENGTH_TIME_RELATION_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3117"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3117_accepted_Omega_dim_sources"], 0)
        self.assertEqual(cert["candidate_R_dim_relations"], 11)
        self.assertEqual(cert["relation_law_rows"], 88)
        self.assertEqual(cert["scale_covariance_rows"], 55)
        self.assertEqual(cert["phase_area_coupling_rows"], 55)
        self.assertEqual(cert["candidate_gate_rows"], 143)
        self.assertEqual(cert["accepted_R_dim_relations"], 0)

    def test_candidate_matrix_and_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_R_dim_relations"]
        self.assertTrue(any(row["candidate"] == "phase_tick_product_relation" and row["C_phi_A_phi_preserved"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "quotient_orbit_section_relation" and not row["not_unit_convention"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "planck_hbar_relation" and not row["standard_physics_import_free"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_oriented_relation" and not row["selector_bridge_ltotal_toe_free"] for row in candidates))
        self.assertTrue(all(not row["accepted_R_dim_relation"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["accepted_coupling_chain"] for row in objs["phase_area_coupling_rows"]))

    def test_negative_flags_docs_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["candidate_R_dim_relations_constructed"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("strict axis-source object Xi_LT", decision["next_honest_step"])
        self.assertIn("distinct length/time axes", decision["next_honest_step"])
        self.assertIn("P3118/S2068", MD.read_text(encoding="utf-8"))
        self.assertIn("P3118/S2068", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3118/S2068", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3118", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
