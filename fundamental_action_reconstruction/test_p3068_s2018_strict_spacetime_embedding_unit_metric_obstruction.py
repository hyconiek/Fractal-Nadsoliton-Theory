import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3068_s2018_strict_spacetime_embedding_unit_metric_obstruction.py"
OUT = ROOT / "generated" / "p3068_s2018_strict_spacetime_embedding_unit_metric_obstruction.json"
MD = ROOT / "generated" / "p3068_s2018_strict_spacetime_embedding_unit_metric_obstruction.md"

class P3068StrictSpacetimeEmbeddingUnitMetricTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input_hash(self):
        self.assertEqual(self.payload["status"], "P3068_STRICT_SPACETIME_EMBEDDING_UNIT_METRIC_OBSTRUCTION_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3067"])

    def test_content_first_and_sat_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["required_embedding_atoms"], 7)
        self.assertEqual(cert["sat_rows"], 128)
        self.assertEqual(cert["accepted_sat_rows"], 1)
        self.assertEqual(cert["minimal_accepted_atom_count"], 7)
        self.assertEqual(cert["provider_candidates"], 5)
        self.assertEqual(cert["accepted_provider_rows"], 0)
        self.assertEqual(cert["current_present_atoms"], 2)
        self.assertEqual(cert["current_missing_atoms"], 5)
        self.assertFalse(cert["current_artifact_accepted"])
        self.assertEqual(cert["p3067_proxy_null_covariance_pass_rows"], 8)
        self.assertEqual(cert["p3067_strict_lorentz_closure_rows"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 4)

    def test_constructed_embedding_template(self):
        obj = self.payload["constructed_theoretical_objects"]["strict_nadsoliton_spacetime_embedding_template"]
        self.assertEqual(obj["object"], "StrictNadsolitonSpacetimeEmbeddingTemplate")
        self.assertIn("E_sigma", obj["formal_shape"])
        self.assertEqual(len(obj["required_atoms"]), 7)
        current = self.payload["constructed_theoretical_objects"]["current_artifact_row"]
        self.assertIn("time_coordinate_map_t", current["missing_atoms"])
        self.assertIn("unit_normalized_metric_g_c", current["missing_atoms"])
        self.assertFalse(current["accepted_strict_embedding"])

    def test_flags_recommendation_and_docs(self):
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["embedding_sat_table_exhausted"])
        self.assertIn("coordinate-pair source theorem", self.payload["decision"]["next_honest_step"])
        self.assertIn("P3068/S2018", MD.read_text(encoding="utf-8"))
        self.assertIn("P3068/S2018", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3068/S2018", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3068", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
