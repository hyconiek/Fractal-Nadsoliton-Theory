import json, subprocess, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3062_s2012_sigma_selector_source_law_candidate_audit.py"
OUT = ROOT / "generated" / "p3062_s2012_sigma_selector_source_law_candidate_audit.json"
MD = ROOT / "generated" / "p3062_s2012_sigma_selector_source_law_candidate_audit.md"

class P3062SigmaSelectorSourceLawCandidateAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_input(self):
        self.assertEqual(self.payload["status"], "P3062_SIGMA_SELECTOR_SOURCE_LAW_CANDIDATE_AUDIT_NO_EXPORT")
        self.assertIsNotNone(self.payload["input_hashes"]["P3061"])

    def test_finite_certificate(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 3)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["acceptance_criteria"], 5)
        self.assertEqual(cert["candidate_classes_audited"], 7)
        self.assertEqual(cert["candidates_with_nonzero_signed_value"], 4)
        self.assertEqual(cert["candidates_coupled_to_G_selector"], 0)
        self.assertEqual(cert["accepted_sigma_selector_source_laws"], 0)
        self.assertEqual(cert["satisfied_proof_obligations"], 3)

    def test_acceptance_boundary_and_candidates(self):
        obj = self.payload["constructed_theoretical_objects"]
        boundary = obj["sigma_selector_source_law_acceptance_boundary"]
        self.assertEqual(boundary["object"], "ConcreteSigmaSelectorSourceLawAcceptanceBoundary")
        self.assertEqual(len(boundary["criteria"]), 5)
        rows = obj["candidate_rows"]
        self.assertEqual(len(rows), 7)
        self.assertTrue(all(not row["accepted_as_sigma_selector_source_law"] for row in rows))
        self.assertTrue(any(row["criteria"]["nonzero_signed_value_computed"] for row in rows))
        for value in self.payload["decision"]["negative_export_flags"].values():
            self.assertFalse(value)

    def test_docs(self):
        self.assertIn("P3062/S2012", MD.read_text(encoding="utf-8"))
        self.assertIn("P3062/S2012", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3062/S2012", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3062", (REPO / "AGENTS.md").read_text(encoding="utf-8"))

if __name__ == "__main__":
    unittest.main()
