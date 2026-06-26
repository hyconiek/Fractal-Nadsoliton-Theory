import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3115_s2065_sigma_dim_scale_section_theorem_audit.py"
OUT = ROOT / "generated" / "p3115_s2065_sigma_dim_scale_section_theorem_audit.json"
MD = ROOT / "generated" / "p3115_s2065_sigma_dim_scale_section_theorem_audit.md"


class P3115SigmaDimScaleSectionTheoremAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3115_SIGMA_DIM_SCALE_SECTION_THEOREM_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3114"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3114_accepted_D_phi_source_laws"], 0)
        self.assertEqual(cert["candidate_Sigma_dim_theorems"], 8)
        self.assertEqual(cert["scale_orbit_witness_rows"], 40)
        self.assertEqual(cert["section_obligation_rows"], 48)
        self.assertEqual(cert["coupling_residual_rows"], 8)
        self.assertEqual(cert["candidate_gate_rows"], 88)
        self.assertEqual(cert["accepted_Sigma_dim_theorems"], 0)

    def test_candidate_matrix_and_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_Sigma_dim_theorems"]
        self.assertTrue(any(row["candidate"] == "alpha_geo_phase_area_section" and row["C_phi_A_phi_coupling_proved"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "unit_norm_triad_section" and row["unique_section_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "planck_light_cone_section" and not row["standard_physics_import_free"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_oriented_dimension_section" and not row["selector_bridge_ltotal_toe_free"] for row in candidates))
        self.assertTrue(all(not row["accepted_Sigma_dim_theorem"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["all_residuals_zero_import_free"] for row in objs["coupling_residual_rows"]))

    def test_negative_flags_docs_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["candidate_Sigma_dim_theorems_constructed"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("dimension-source functor K_dim", decision["next_honest_step"])
        self.assertIn("positive scale torsor", decision["next_honest_step"])
        self.assertIn("P3115/S2065", MD.read_text(encoding="utf-8"))
        self.assertIn("P3115/S2065", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3115/S2065", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3115", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
