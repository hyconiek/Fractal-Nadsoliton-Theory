import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3150_s2100_rsm_hypercharge_ray_source_selection_audit.py"
OUT = ROOT / "generated" / "p3150_s2100_rsm_hypercharge_ray_source_selection_audit.json"
MD = ROOT / "generated" / "p3150_s2100_rsm_hypercharge_ray_source_selection_audit.md"


class P3150RsmHyperchargeRaySourceSelectionAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_linear_counts(self):
        self.assertEqual(self.payload["status"], "P3150_RSM_HYPERCHARGE_RAY_SOURCE_SELECTION_CONDITIONAL_RAY_NO_STRICT_SOURCE")
        counts = self.payload["finite_theorem"]["finite_counts"]
        self.assertEqual(counts["linear_equation_rows"], 5)
        self.assertEqual(counts["unknowns"], 6)
        self.assertEqual(counts["matrix_rank"], 5)
        self.assertEqual(counts["nullity"], 1)
        self.assertEqual(counts["ray_witnesses_matching_P3148"], 1)
        self.assertEqual(counts["redundant_anomaly_checks_vanishing"], 1)
        self.assertEqual(counts["strict_accepted_source_rows"], 0)

    def test_normalized_ray_matches_standard_assignments(self):
        ray = self.payload["hypercharge_ray_derivation"]
        self.assertTrue(ray["matches_P3148_standard_hypercharges"])
        self.assertTrue(ray["redundant_checks_vanish"])
        self.assertEqual(ray["normalized_h_equals_one_half"], {
            "q": "1/6",
            "u": "-2/3",
            "d": "1/3",
            "l": "-1/2",
            "e": "1",
            "h": "1/2",
        })

    def test_candidate_rows_do_not_export_strict_source(self):
        candidates = self.payload["candidate_source_rows"]
        self.assertEqual(len(candidates), 4)
        self.assertTrue(any(row["candidate"] == "Y_SM^ray consistency constraints" and row["selects_hypercharge_ray"] for row in candidates))
        self.assertFalse(any(row["strict_nadsoliton_source_law"] and row["noncircular_current_artifact"] for row in candidates))
        decision = self.payload["decision"]
        self.assertIn("P3151", decision["next_honest_step"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))

    def test_docs_updated(self):
        self.assertIn("P3150/S2100", MD.read_text(encoding="utf-8"))
        self.assertIn("P3150/S2100", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3150/S2100", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3150", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
