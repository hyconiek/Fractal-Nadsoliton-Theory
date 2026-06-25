import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3084_s2034_gauge_representation_obstruction_witness_audit.py"
OUT = ROOT / "generated" / "p3084_s2034_gauge_representation_obstruction_witness_audit.json"
MD = ROOT / "generated" / "p3084_s2034_gauge_representation_obstruction_witness_audit.md"

class P3084GaugeRepresentationObstructionWitnessAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3084_GAUGE_REPRESENTATION_OBSTRUCTION_WITNESS_BOUNDED_NO_GO")
        self.assertIsNotNone(self.payload["input_hashes"]["P3083"])

    def test_finite_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3083_accepted_nonimported_lorentz_signature_sources"], 0)
        self.assertEqual(cert["z12_character_rows"], 12)
        self.assertEqual(cert["nontrivial_character_rows"], 11)
        self.assertEqual(cert["character_orthogonality_rows"], 144)
        self.assertEqual(cert["orthogonality_failures"], 0)
        self.assertEqual(cert["flat_holonomy_rows"], 12)
        self.assertEqual(cert["nonzero_curvature_rows"], 0)
        self.assertEqual(cert["twisted_laplacian_rows"], 12)
        self.assertEqual(cert["twisted_rows_with_sourced_gauge_dynamics"], 0)
        self.assertEqual(cert["gauge_candidates"], 5)
        self.assertEqual(cert["candidate_gate_rows"], 30)
        self.assertEqual(cert["accepted_nonimported_gauge_representation_sources"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_character_and_candidate_aggregates(self):
        chars = self.payload["constructed_theoretical_objects"]["z12_character_rows"]
        self.assertEqual(chars[0]["order_divisor"], 1)
        self.assertFalse(chars[0]["nontrivial"])
        self.assertTrue(chars[1]["nontrivial"])
        aggs = self.payload["constructed_theoretical_objects"]["candidate_aggregate_certificate"]
        self.assertTrue(all(not row["accepted_nonimported_gauge_representation_source"] for row in aggs))
        by_id = {row["candidate"]: row for row in aggs}
        self.assertEqual(by_id["z12_regular_character_table"]["passed_gates"], 1)
        self.assertEqual(by_id["flat_cycle_holonomy_labels"]["passed_gates"], 2)
        self.assertEqual(by_id["imported_u1_minimal_coupling_template"]["passed_gates"], 5)

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["character_orthogonality_verified"])
        self.assertIn("conserved-current/Noether-obstruction", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3084/S2034", MD.read_text(encoding="utf-8"))
        self.assertIn("P3084/S2034", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3084/S2034", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3084", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
