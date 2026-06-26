import json, subprocess, sys, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3109_s2059_alpha_geo_scale_orbit_quotient_source_law_audit.py"
OUT = ROOT / "generated" / "p3109_s2059_alpha_geo_scale_orbit_quotient_source_law_audit.json"
MD = ROOT / "generated" / "p3109_s2059_alpha_geo_scale_orbit_quotient_source_law_audit.md"

class P3109AlphaGeoScaleOrbitQuotientSourceLawAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3109_ALPHA_GEO_SCALE_ORBIT_QUOTIENT_SOURCE_LAW_OBSTRUCTION_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3108"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3108_accepted_internal_calibration_morphisms"], 0)
        self.assertEqual(cert["entropy_cell_rows"], 4)
        self.assertEqual(cert["alpha_geo_four_bit_cells"], 1)
        self.assertEqual(cert["positive_scale_orbit_rows"], 9)
        self.assertEqual(cert["section_candidates"], 5)
        self.assertEqual(cert["calibration_coupling_rows"], 15)
        self.assertEqual(cert["candidate_gate_rows"], 40)
        self.assertEqual(cert["accepted_scale_orbit_quotient_source_laws"], 0)

    def test_quotient_class_positive_but_no_section(self):
        objs = self.payload["constructed_theoretical_objects"]
        self.assertIn("self-coupled primordial information", objs["quotient_source_law_audit_object"]["ontology"])
        self.assertTrue(any(row["is_alpha_geo_cell"] and row["states"] == 16 for row in objs["entropy_cell_rows"]))
        self.assertTrue(all(row["quotient_class_fixed_by_entropy"] for row in objs["positive_scale_orbit_rows"]))
        self.assertTrue(all(not row["unique_section_selected"] for row in objs["positive_scale_orbit_rows"]))
        self.assertTrue(all(not row["accepted_scale_orbit_quotient_source_law"] for row in objs["candidate_aggregate_certificate"]))

    def test_negative_flags_docs_and_recommendation(self):
        decision = self.payload["decision"]
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("dimension-carrying comparison standard", decision["next_honest_step"])
        self.assertIn("P3109/S2059", MD.read_text(encoding="utf-8"))
        self.assertIn("P3109/S2059", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3109/S2059", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3109", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
