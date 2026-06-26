import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3116_s2066_k_dim_dimension_source_functor_audit.py"
OUT = ROOT / "generated" / "p3116_s2066_k_dim_dimension_source_functor_audit.json"
MD = ROOT / "generated" / "p3116_s2066_k_dim_dimension_source_functor_audit.md"


class P3116KDimDimensionSourceFunctorAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3116_K_DIM_DIMENSION_SOURCE_FUNCTOR_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3115"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3115_accepted_Sigma_dim_theorems"], 0)
        self.assertEqual(cert["candidate_K_dim_functors"], 9)
        self.assertEqual(cert["functor_law_rows"], 54)
        self.assertEqual(cert["torsor_action_rows"], 45)
        self.assertEqual(cert["source_coupling_residual_rows"], 36)
        self.assertEqual(cert["candidate_gate_rows"], 99)
        self.assertEqual(cert["accepted_K_dim_functors"], 0)

    def test_candidate_matrix_and_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_K_dim_functors"]
        self.assertTrue(any(row["candidate"] == "alpha_geo_phase_area_functor" and row["C_phi_A_phi_coupling_sourced"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "constant_unit_torsor_functor" and row["functor_laws_verified"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "planck_unit_functor" and not row["standard_physics_import_free"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_choice_functor" and not row["selector_bridge_ltotal_toe_free"] for row in candidates))
        self.assertTrue(all(not row["accepted_K_dim_functor"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["accepted_zero"] for row in objs["source_coupling_residual_rows"]))

    def test_negative_flags_docs_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["candidate_K_dim_functors_constructed"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("strict typed source object Omega_dim", decision["next_honest_step"])
        self.assertIn("internal dimension character/valuation", decision["next_honest_step"])
        self.assertIn("P3116/S2066", MD.read_text(encoding="utf-8"))
        self.assertIn("P3116/S2066", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3116/S2066", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3116", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
