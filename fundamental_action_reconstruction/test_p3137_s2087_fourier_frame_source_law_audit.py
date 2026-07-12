import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3137_s2087_fourier_frame_source_law_audit.py"
OUT = ROOT / "generated" / "p3137_s2087_fourier_frame_source_law_audit.json"
MD = ROOT / "generated" / "p3137_s2087_fourier_frame_source_law_audit.md"


class P3137FourierFrameSourceLawAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_backscan(self):
        self.assertEqual(self.payload["status"], "P3137_FOURIER_FRAME_SOURCE_LAW_F_DHL_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3136"])
        self.assertEqual(len(self.payload["repo_backscan_summary"]), 4)

    def test_finite_frame_certificate(self):
        cert = self.payload["finite_frame_certificate"]
        self.assertEqual(cert["profiles_tested"], 120)
        self.assertEqual(cert["primitive_character_orbit_size"], 4)
        self.assertEqual(cert["primitive_character_orbit"], [1, 5, 7, 11])
        self.assertEqual(cert["primitive_pair_orbit_size"], 2)
        self.assertEqual(cert["profiles_with_primitive_active_pair"], 48)
        self.assertEqual(cert["profiles_with_nonprimitive_active_pair"], 72)
        self.assertEqual(cert["phase_zero_cells"], 12)
        self.assertEqual(cert["source_candidates_tested"], 9)
        self.assertEqual(cert["source_candidates_passing_all_gates"], 0)
        self.assertEqual(cert["accepted_F_DHL_sources"], 0)

    def test_source_rows(self):
        rows = self.payload["source_acceptance_rows"]
        self.assertTrue(any(row["candidate_F_DHL_source"] == "active_spectral_pair_receiver" and row["selects_character_or_pair"] for row in rows))
        self.assertTrue(any(row["candidate_F_DHL_source"] == "P2992_frequency_source_localizer" and row["import_free_current_artifacts"] for row in rows))
        self.assertTrue(any(row["candidate_F_DHL_source"] == "phase_zero_by_argument_gauge" and row["selects_phase_zero_reference"] and not row["import_free_current_artifacts"] for row in rows))
        self.assertTrue(all(not row["accepted_F_DHL_source"] for row in rows))

    def test_decision_and_docs(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["F_DHL_candidate_space_constructed"])
        self.assertTrue(decision["positive_scoped_flags"]["primitive_character_orbit_obstruction_proved"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("non-Fourier joint source object", decision["next_honest_step"])
        self.assertIn("P3137/S2087", MD.read_text(encoding="utf-8"))
        self.assertIn("P3137/S2087", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3137/S2087", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3137", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
