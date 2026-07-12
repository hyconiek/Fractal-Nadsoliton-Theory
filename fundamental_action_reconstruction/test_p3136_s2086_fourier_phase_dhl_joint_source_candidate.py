import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3136_s2086_fourier_phase_dhl_joint_source_candidate.py"
OUT = ROOT / "generated" / "p3136_s2086_fourier_phase_dhl_joint_source_candidate.json"
MD = ROOT / "generated" / "p3136_s2086_fourier_phase_dhl_joint_source_candidate.md"


class P3136FourierPhaseJDHLCandidateAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_constructed_formula(self):
        self.assertEqual(self.payload["status"], "P3136_FOURIER_PHASE_J_DHL_CANDIDATE_CHART_RELATIVE_POSITIVE_STRICT_SOURCE_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3135"])
        self.assertIn("C_k(D)", self.payload["constructed_object"]["formula"])
        self.assertGreaterEqual(len(self.payload["repo_backscan_summary"]), 3)

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["candidate_profiles_tested"], 120)
        self.assertEqual(cert["nonzero_fourier_coefficients"], 120)
        self.assertEqual(cert["unit_mode_profiles"], 48)
        self.assertEqual(cert["nonunit_mode_profiles"], 72)
        self.assertEqual(cert["chart_relative_full_pair_recoveries"], 0)
        self.assertEqual(cert["aliased_nonunit_recoveries"], 120)
        self.assertEqual(cert["translation_phase_rotation_witnesses"], 120)
        self.assertEqual(cert["accepted_import_free_J_DHL_sources"], 0)
        self.assertEqual(cert["mode_alias_counts"], {"1": [2], "2": [4], "3": [6], "4": [4], "5": [2]})

    def test_rows_capture_aliasing_and_unit_recovery(self):
        reps = {row["mode_k"]: row for row in self.payload["representative_rows"]}
        self.assertFalse(reps[1]["chart_relative_full_pair_recovery"])
        self.assertEqual(reps[1]["matching_pairs"], [[0, 1], [6, -1]])
        self.assertFalse(reps[2]["chart_relative_full_pair_recovery"])
        self.assertEqual(len(reps[2]["matching_pairs"]), 4)
        self.assertEqual(len(reps[3]["matching_pairs"]), 6)
        self.assertEqual(len(reps[4]["matching_pairs"]), 4)
        self.assertFalse(reps[5]["chart_relative_full_pair_recovery"])
        self.assertTrue(all(row["translation_phase_rotates_at_t1"] for row in self.payload["extraction_rows"]))
        self.assertTrue(all(not row["strict_import_free_recovery"] for row in self.payload["extraction_rows"]))

    def test_decision_and_docs(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["formula_level_J_DHL_candidate_constructed"])
        self.assertTrue(decision["positive_scoped_flags"]["unit_mode_half_period_aliasing_verified"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("F_DHL", decision["next_honest_step"])
        self.assertIn("P3136/S2086", MD.read_text(encoding="utf-8"))
        self.assertIn("P3136/S2086", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3136/S2086", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3136", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
