import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3163_s2113_boundary_value_speed_of_light_matching_audit.py"
OUT = ROOT / "generated" / "p3163_s2113_boundary_value_speed_of_light_matching_audit.json"
MD = ROOT / "generated" / "p3163_s2113_boundary_value_speed_of_light_matching_audit.md"


class P3163BoundaryValueSpeedOfLightMatchingAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_boundary_fit_counts(self):
        self.assertEqual(self.payload["status"], "P3163_BOUNDARY_VALUE_SPEED_OF_LIGHT_MATCHING_UNDERDETERMINATION_AUDIT")
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["boundary_value_candidates"], 10)
        self.assertEqual(cert["speed_of_light_fit_rows"], 10)
        self.assertEqual(cert["numerically_fit_nonzero_candidates"], 9)
        self.assertEqual(cert["scale_degeneracy_rows"], 225)
        self.assertEqual(cert["accepted_physical_c_sources"], 0)
        self.assertEqual(cert["c_SI_m_per_s"], 299792458.0)

    def test_fit_rows_are_not_sources(self):
        objs = self.payload["constructed_theoretical_objects"]
        fits = objs["speed_of_light_fit_rows"]
        self.assertTrue(any(row["candidate"] == "strict_kernel_tail_limit" and not row["can_numerically_fit_c"] for row in fits))
        self.assertTrue(all(not row["source_closure"] for row in fits))
        self.assertTrue(all(not row["accepted_physical_c_source"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(value is False for value in self.payload["decision"]["negative_export_flags"].values()))
        self.assertIn("U_length", self.payload["decision"]["next_honest_step"])
        self.assertIn("U_time", self.payload["decision"]["next_honest_step"])

    def test_docs_updated(self):
        self.assertIn("P3163/S2113", MD.read_text(encoding="utf-8"))
        self.assertIn("P3163/S2113", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3163/S2113", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3163", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
