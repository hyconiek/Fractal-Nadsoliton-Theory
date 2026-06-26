import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3117_s2067_omega_dim_dimension_character_source_audit.py"
OUT = ROOT / "generated" / "p3117_s2067_omega_dim_dimension_character_source_audit.json"
MD = ROOT / "generated" / "p3117_s2067_omega_dim_dimension_character_source_audit.md"


class P3117OmegaDimDimensionCharacterSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3117_OMEGA_DIM_DIMENSION_CHARACTER_SOURCE_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3116"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3116_accepted_K_dim_functors"], 0)
        self.assertEqual(cert["candidate_Omega_dim_characters"], 10)
        self.assertEqual(cert["scale_valuation_rows"], 60)
        self.assertEqual(cert["dimension_axis_rows"], 50)
        self.assertEqual(cert["coupling_chain_rows"], 40)
        self.assertEqual(cert["candidate_gate_rows"], 120)
        self.assertEqual(cert["accepted_Omega_dim_sources"], 0)

    def test_candidate_matrix_and_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_Omega_dim_characters"]
        self.assertTrue(any(row["candidate"] == "phase_area_character" and row["sources_C_phi_A_phi_equals_U_action"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "internal_tick_ratio_character" and not row["not_gauge_convention"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "planck_dimensional_character" and not row["standard_physics_import_free"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_oriented_character" and not row["selector_bridge_ltotal_toe_free"] for row in candidates))
        self.assertTrue(all(not row["accepted_Omega_dim_source"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["accepted_coupling"] for row in objs["coupling_chain_rows"]))

    def test_negative_flags_docs_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["candidate_Omega_dim_characters_constructed"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("strict relation object R_dim", decision["next_honest_step"])
        self.assertIn("action-length-time composition law", decision["next_honest_step"])
        self.assertIn("P3117/S2067", MD.read_text(encoding="utf-8"))
        self.assertIn("P3117/S2067", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3117/S2067", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3117", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
