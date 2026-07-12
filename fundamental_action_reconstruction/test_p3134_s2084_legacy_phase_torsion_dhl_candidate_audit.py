import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3134_s2084_legacy_phase_torsion_dhl_candidate_audit.py"
OUT = ROOT / "generated" / "p3134_s2084_legacy_phase_torsion_dhl_candidate_audit.json"
MD = ROOT / "generated" / "p3134_s2084_legacy_phase_torsion_dhl_candidate_audit.md"


class P3134LegacyPhaseTorsionDHLCandidateAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_inputs_and_formula(self):
        self.assertEqual(self.payload["status"], "P3134_EXPLICIT_LEGACY_PHASE_TORSION_DHL_CANDIDATE_ORIGIN_POLARITY_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3133"])
        self.assertIn("sin(2*pi*k*(x-r)/12)", self.payload["constructed_object"]["formula"])

    def test_finite_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["origins_r"], 12)
        self.assertEqual(cert["nonzero_modes_k"], 5)
        self.assertEqual(cert["lambda_polarities"], 2)
        self.assertEqual(cert["candidate_defects"], 120)
        self.assertEqual(cert["nonzero_odd_candidates"], 120)
        self.assertEqual(cert["support_coupled_if_origin_given"], 120)
        self.assertEqual(cert["translation_symmetry_rows_sampled"], 240)
        self.assertEqual(cert["equivariant_with_origin_shift_rows"], 240)
        self.assertEqual(cert["quotient_invariant_nonzero_t1_rows"], 0)
        self.assertEqual(cert["accepted_import_free_D_HL_sources"], 0)

    def test_rows_show_conditional_object_not_export(self):
        rows = self.payload["candidate_rows"]
        self.assertTrue(all(row["nonzero_defect"] for row in rows))
        self.assertTrue(all(row["odd_around_chosen_r"] for row in rows))
        self.assertTrue(all(row["support_coupled_if_r_is_given"] for row in rows))
        self.assertTrue(all(not row["accepted_import_free_D_HL"] for row in rows))
        t1 = [row for row in self.payload["symmetry_rows_sampled"] if row["translation_t"] == 1]
        self.assertEqual(len(t1), 120)
        self.assertTrue(all(row["equivariant_if_origin_shifts"] for row in t1))
        self.assertTrue(all(not row["invariant_if_origin_quotiented"] for row in t1))

    def test_decision_and_docs(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["explicit_D_HL_formula_constructed"])
        self.assertTrue(decision["positive_scoped_flags"]["origin_and_lambda_missing_premises_isolated"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("joint (r,lambda) selector-source matrix", decision["next_honest_step"])
        self.assertIn("P3134/S2084", MD.read_text(encoding="utf-8"))
        self.assertIn("P3134/S2084", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3134/S2084", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3134", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
