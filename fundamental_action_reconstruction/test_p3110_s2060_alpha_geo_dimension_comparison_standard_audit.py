import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3110_s2060_alpha_geo_dimension_comparison_standard_audit.py"
OUT = ROOT / "generated" / "p3110_s2060_alpha_geo_dimension_comparison_standard_audit.json"
MD = ROOT / "generated" / "p3110_s2060_alpha_geo_dimension_comparison_standard_audit.md"

class P3110AlphaGeoDimensionComparisonStandardAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3110_ALPHA_GEO_DIMENSION_COMPARISON_STANDARD_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3109"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3109_accepted_scale_orbit_quotient_source_laws"], 0)
        self.assertEqual(cert["candidate_standards"], 6)
        self.assertEqual(cert["targets"], 3)
        self.assertEqual(cert["scale_orbit_section_rows"], 18)
        self.assertEqual(cert["required_gates"], 8)
        self.assertEqual(cert["candidate_gate_rows"], 48)
        self.assertEqual(cert["accepted_dimension_comparison_standards"], 0)

    def test_constructed_objects_and_best_candidate_boundary(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertIn("self-coupled primordial information", objs["audit_object"]["ontology"])
        self.assertTrue(any(row["candidate"] == "symplectic_phase_area_candidate" and row["dimension_bearing_standard"] for row in objs["candidate_standard_rows"]))
        self.assertTrue(all(not row["accepted_dimension_comparison_standard"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["unique_positive_section_selected"] for row in objs["scale_orbit_section_test_rows"] if row["candidate"] != "imported_planck_action_comparison"))

    def test_negative_flags_docs_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("symplectic/action phase source-law", decision["next_honest_step"])
        self.assertIn("P3110/S2060", MD.read_text(encoding="utf-8"))
        self.assertIn("P3110/S2060", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3110/S2060", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3110", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
