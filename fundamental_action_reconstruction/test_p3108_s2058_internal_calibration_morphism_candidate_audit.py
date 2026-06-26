import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3108_s2058_internal_calibration_morphism_candidate_audit.py"
OUT = ROOT / "generated" / "p3108_s2058_internal_calibration_morphism_candidate_audit.json"
MD = ROOT / "generated" / "p3108_s2058_internal_calibration_morphism_candidate_audit.md"

class P3108InternalCalibrationMorphismCandidateAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3108_INTERNAL_CALIBRATION_MORPHISM_CANDIDATE_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3107"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3107_accepted_internal_shannon_to_dimension_source_laws"], 0)
        self.assertEqual(cert["calibration_objects"], 4)
        self.assertEqual(cert["morphism_candidates"], 6)
        self.assertEqual(cert["scale_orbit_rows"], 24)
        self.assertEqual(cert["source_law_rows"], 6)
        self.assertEqual(cert["candidate_gate_rows"], 54)
        self.assertEqual(cert["accepted_internal_calibration_morphisms"], 0)

    def test_positive_internal_entropy_but_no_units(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertIn("self-coupled primordial information", objs["calibration_morphism_audit_object"]["ontology"])
        self.assertTrue(any(row["source_object"] == "I4_alpha_geo_shannon_cell" for row in objs["morphism_candidate_rows"]))
        self.assertTrue(any(row["internal_formula_defined"] for row in objs["morphism_candidate_rows"]))
        self.assertTrue(all(not row["passes_scale_orbit_test"] for row in objs["scale_orbit_rows"] if row["target_nonzero_dimensional"] or row["scale_breaking_claimed"]))
        self.assertTrue(all(not row["accepted_internal_calibration_morphism"] for row in objs["candidate_aggregate_certificate"]))

    def test_negative_flags_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("scale-orbit quotient", decision["next_honest_step"])
        self.assertIn("P3108/S2058", MD.read_text(encoding="utf-8"))
        self.assertIn("P3108/S2058", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3108/S2058", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3108", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
