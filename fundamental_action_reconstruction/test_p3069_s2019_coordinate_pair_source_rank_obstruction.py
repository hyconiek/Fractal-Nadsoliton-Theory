import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3069_s2019_coordinate_pair_source_rank_obstruction.py"
OUT = ROOT / "generated" / "p3069_s2019_coordinate_pair_source_rank_obstruction.json"
MD = ROOT / "generated" / "p3069_s2019_coordinate_pair_source_rank_obstruction.md"

class P3069CoordinatePairSourceRankObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3069_COORDINATE_PAIR_SOURCE_RANK_OBSTRUCTION_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3068"])

    def test_content_first_and_rank_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["domain_rows"], 24)
        self.assertEqual(cert["candidate_features"], 6)
        self.assertEqual(cert["feature_subset_rows"], 63)
        self.assertGreater(cert["raw_rank_two_or_more_rows"], 0)
        self.assertEqual(cert["accepted_coordinate_pair_source_rows"], 0)
        self.assertEqual(cert["unit_bearing_candidate_features"], 0)
        self.assertEqual(cert["nonconventional_source_features"], 0)
        self.assertGreaterEqual(cert["max_raw_rank"], 2)
        self.assertEqual(cert["p3068_current_missing_atoms"], 5)
        self.assertFalse(cert["p3068_current_artifact_accepted"])
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_constructed_candidate_space(self):
        obj = self.payload["constructed_theoretical_objects"]["coordinate_pair_source_candidate_space"]
        self.assertEqual(obj["object"], "CoordinatePairSourceCandidateSpace")
        self.assertEqual(len(obj["candidate_features"]), 6)
        matrix = self.payload["constructed_theoretical_objects"]["coordinate_pair_rank_obstruction_matrix"]
        self.assertEqual(len(matrix), 63)
        self.assertTrue(any(row["raw_two_column_span"] for row in matrix))
        self.assertTrue(all(not row["accepted_coordinate_pair_source"] for row in matrix))

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["rational_rank_matrix_executed"])
        self.assertIn("canonical length/time unit provider", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3069/S2019", MD.read_text(encoding="utf-8"))
        self.assertIn("P3069/S2019", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3069/S2019", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3069", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
