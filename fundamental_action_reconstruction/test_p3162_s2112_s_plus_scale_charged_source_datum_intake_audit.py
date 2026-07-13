import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3162_s2112_s_plus_scale_charged_source_datum_intake_audit.py"
OUT = ROOT / "generated" / "p3162_s2112_s_plus_scale_charged_source_datum_intake_audit.json"
MD = ROOT / "generated" / "p3162_s2112_s_plus_scale_charged_source_datum_intake_audit.md"


class P3162SPlusScaleChargedSourceDatumIntakeAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_s_plus_object_and_counts(self):
        self.assertEqual(self.payload["status"], "P3162_S_PLUS_SCALE_CHARGED_SOURCE_DATUM_INTAKE_BOUNDED_NO_GO")
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["candidate_S_plus_sources"], 12)
        self.assertEqual(cert["scale_factors"], 5)
        self.assertEqual(cert["representation_rows"], 60)
        self.assertEqual(cert["candidate_gate_rows"], 108)
        self.assertEqual(cert["accepted_S_plus_sources"], 0)
        self.assertEqual(cert["accepted_weight_rows"], 0)
        self.assertGreater(cert["weight_plus_one_candidates"], 0)
        self.assertGreater(cert["weight_zero_candidates"], 0)

    def test_weight_one_schema_but_no_export(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_S_plus_sources"]
        self.assertTrue(any(row["candidate"] == "new_strict_scale_charged_datum_schema" and row["scale_weight"] == 1 for row in candidates))
        self.assertTrue(any(row["candidate"] == "alpha_geo_phase_area" and row["scale_weight"] == 0 for row in candidates))
        self.assertTrue(all(not row["accepted_S_plus_source"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))
        self.assertTrue(self.payload["decision"]["positive_scoped_flags"]["weight_one_acceptance_schema_defined"])

    def test_docs_updated(self):
        self.assertIn("P3162/S2112", MD.read_text(encoding="utf-8"))
        self.assertIn("P3162/S2112", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3162/S2112", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3162", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
