import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3070_s2020_canonical_length_time_unit_provider_obstruction.py"
OUT = ROOT / "generated" / "p3070_s2020_canonical_length_time_unit_provider_obstruction.json"
MD = ROOT / "generated" / "p3070_s2020_canonical_length_time_unit_provider_obstruction.md"

class P3070CanonicalLengthTimeUnitProviderObstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3070_CANONICAL_LENGTH_TIME_UNIT_PROVIDER_OBSTRUCTION_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3069"])

    def test_content_first_and_matrix_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["candidate_scalars"], 6)
        self.assertEqual(cert["candidate_unit_maps"], 5)
        self.assertEqual(cert["provider_matrix_rows"], 30)
        self.assertGreater(cert["positive_scalar_candidates"], 0)
        self.assertGreater(cert["intrinsic_scalar_candidates"], 0)
        self.assertEqual(cert["nonconventional_scalar_candidates"], 0)
        self.assertEqual(cert["unit_maps_with_length_and_time"], 1)
        self.assertEqual(cert["unit_maps_with_coordinate_coupling"], 1)
        self.assertEqual(cert["nonconventional_unit_laws"], 0)
        self.assertEqual(cert["accepted_unit_provider_rows"], 0)
        self.assertEqual(cert["p3069_accepted_coordinate_pair_source_rows"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_constructed_template_and_obstruction_rows(self):
        obj = self.payload["constructed_theoretical_objects"]["canonical_unit_provider_template"]
        self.assertEqual(obj["object"], "CanonicalLengthTimeUnitProviderTemplate")
        self.assertIn("ell_s", obj["formal_shape"])
        rows = self.payload["constructed_theoretical_objects"]["unit_provider_obstruction_matrix"]
        self.assertEqual(len(rows), 30)
        self.assertTrue(any(row["positive_scale"] for row in rows))
        self.assertTrue(all(not row["accepted_canonical_length_time_unit_provider"] for row in rows))

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["finite_scalar_unit_matrix_executed"])
        self.assertIn("sigma-invariant scalar conservation/scale-control", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3070/S2020", MD.read_text(encoding="utf-8"))
        self.assertIn("P3070/S2020", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3070/S2020", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3070", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
