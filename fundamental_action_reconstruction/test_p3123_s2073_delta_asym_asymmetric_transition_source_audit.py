import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3123_s2073_delta_asym_asymmetric_transition_source_audit.py"
OUT = ROOT / "generated" / "p3123_s2073_delta_asym_asymmetric_transition_source_audit.json"
MD = ROOT / "generated" / "p3123_s2073_delta_asym_asymmetric_transition_source_audit.md"


class P3123DeltaAsymAsymmetricTransitionSourceAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3123_DELTA_ASYM_ASYMMETRIC_TRANSITION_SOURCE_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3122"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3122_accepted_Iota_irrev_sources"], 0)
        self.assertEqual(cert["candidate_Delta_asym_sources"], 16)
        self.assertEqual(cert["asymmetry_law_rows"], 192)
        self.assertEqual(cert["transition_witness_rows"], 128)
        self.assertEqual(cert["Iota_Kappa_Tau_Xi_R_coupling_rows"], 144)
        self.assertEqual(cert["candidate_gate_rows"], 288)
        self.assertGreaterEqual(cert["phase_information_coupling_rows"], 1)
        self.assertEqual(cert["accepted_Delta_asym_sources"], 0)

    def test_candidate_matrix_and_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        candidates = objs["candidate_Delta_asym_sources"]
        self.assertTrue(any(row["candidate"] == "phase_information_gradient_asymmetry" and row["nonzero_asymmetry_witness_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "phase_information_gradient_asymmetry" and not row["orientation_not_selector_premise"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "a_phi_weighted_transition_area" and row["C_phi_A_phi_preserved"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "thermodynamic_phase_entropy_production" and not row["strict_nadsoliton_data_only"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_oriented_phase_asymmetry" and not row["Kappa_cycle_source_unlocked"] for row in candidates))
        self.assertTrue(all(not row["accepted_Delta_asym_source"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["accepted_coupling_chain"] for row in objs["Iota_Kappa_Tau_Xi_R_coupling_rows"]))

    def test_phase_information_recommendation_and_docs(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["phase_information_coupling_is_promising_but_scoped"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("information works well with phase", decision["bounded_result"])
        self.assertIn("strict phase-information gauge quotient object Phi_Info", decision["next_honest_step"])
        self.assertIn("phase-origin gauge", decision["next_honest_step"])
        self.assertIn("P3123/S2073", MD.read_text(encoding="utf-8"))
        self.assertIn("P3123/S2073", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3123/S2073", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3123", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
