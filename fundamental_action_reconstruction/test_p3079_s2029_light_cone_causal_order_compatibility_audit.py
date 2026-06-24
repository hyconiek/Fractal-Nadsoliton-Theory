import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3079_s2029_light_cone_causal_order_compatibility_audit.py"
OUT = ROOT / "generated" / "p3079_s2029_light_cone_causal_order_compatibility_audit.json"
MD = ROOT / "generated" / "p3079_s2029_light_cone_causal_order_compatibility_audit.md"

class P3079LightConeCausalOrderCompatibilityAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3079_LIGHT_CONE_CAUSAL_ORDER_COMPATIBILITY_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3078"])

    def test_matrix_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3078_accepted_intrinsic_momentum_symplectic_sources"], 0)
        self.assertEqual(cert["cycle_distance_rows"], 144)
        self.assertEqual(cert["distance_shell_rows"], 7)
        self.assertEqual(cert["distance_shell_counts"], {"0": 12, "1": 24, "2": 24, "3": 24, "4": 24, "5": 24, "6": 12})
        self.assertEqual(cert["diffusion_support_rows"], 21)
        self.assertEqual(cert["diffusion_outside_unit_cone_violations_if_read_as_light"], 15)
        self.assertEqual(cert["spectral_velocity_proxy_rows"], 12)
        self.assertEqual(cert["unit_light_speed_certified_rows"], 0)
        self.assertEqual(cert["causal_candidates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 25)
        self.assertEqual(cert["accepted_internal_causal_order_sources"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 5)

    def test_candidate_aggregates_and_diffusion_rows(self):
        aggs = self.payload["constructed_theoretical_objects"]["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_internal_causal_order_source"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["adjacency_iterate_combinatorial_cone"]["passed_gates"], 2)
        self.assertEqual(by_id["imported_minkowski_light_cone_ansatz"]["passed_gates"], 4)
        diffusion = self.payload["constructed_theoretical_objects"]["diffusion_support_rows"]
        outside = [row for row in diffusion if row["violates_finite_unit_cone_if_read_as_light"]]
        self.assertEqual(len(outside), 15)
        self.assertTrue(all(row["nonzero_tail"] for row in outside))

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["diffusion_tail_obstruction_computed"])
        self.assertIn("typed observable-interface obligation table", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3079/S2029", MD.read_text(encoding="utf-8"))
        self.assertIn("P3079/S2029", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3079/S2029", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3079", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
