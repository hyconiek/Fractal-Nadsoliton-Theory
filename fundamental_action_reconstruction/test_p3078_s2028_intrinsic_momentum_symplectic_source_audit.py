import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3078_s2028_intrinsic_momentum_symplectic_source_audit.py"
OUT = ROOT / "generated" / "p3078_s2028_intrinsic_momentum_symplectic_source_audit.json"
MD = ROOT / "generated" / "p3078_s2028_intrinsic_momentum_symplectic_source_audit.md"

class P3078IntrinsicMomentumSymplecticSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3078_INTRINSIC_MOMENTUM_SYMPLECTIC_SOURCE_NOT_EXPORTED_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3077"])

    def test_matrix_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3077_accepted_internal_second_order_wave_rows"], 0)
        self.assertEqual(cert["candidate_sources"], 5)
        self.assertEqual(cert["required_gates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 25)
        self.assertEqual(cert["z12_skew_current_modal_rows"], 12)
        self.assertEqual(cert["z12_skew_current_zero_rows"], 2)
        self.assertEqual(cert["accepted_intrinsic_momentum_symplectic_sources"], 0)
        self.assertEqual(cert["formal_imported_candidate_passed_gates"], 3)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_candidate_aggregates_and_skew_rows(self):
        aggs = self.payload["constructed_theoretical_objects"]["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_intrinsic_momentum_symplectic_source"] for row in aggs))
        by_id = {row["candidate_source"]: row for row in aggs}
        self.assertEqual(by_id["formal_imported_canonical_phase_space"]["passed_gates"], 3)
        self.assertEqual(by_id["z12_shift_derivative_current"]["passed_gates"], 1)
        skew = self.payload["constructed_theoretical_objects"]["z12_skew_current_modal_rows"]
        self.assertEqual([row["mode_j"] for row in skew if row["zero_skew_eigenvalue"]], [0, 6])
        self.assertTrue(all(not row["can_pair_as_canonical_momentum"] for row in skew))

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["second_order_wave_promotion_frozen_on_current_artifacts"])
        self.assertIn("light-cone/causal-order compatibility audit", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3078/S2028", MD.read_text(encoding="utf-8"))
        self.assertIn("P3078/S2028", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3078/S2028", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3078", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
