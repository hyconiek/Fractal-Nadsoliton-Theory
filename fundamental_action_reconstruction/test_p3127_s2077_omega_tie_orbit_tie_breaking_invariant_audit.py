import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3127_s2077_omega_tie_orbit_tie_breaking_invariant_audit.py"
OUT = ROOT / "generated" / "p3127_s2077_omega_tie_orbit_tie_breaking_invariant_audit.json"
MD = ROOT / "generated" / "p3127_s2077_omega_tie_orbit_tie_breaking_invariant_audit.md"


class P3127OmegaTieOrbitTieBreakingInvariantAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3127_OMEGA_TIE_ORBIT_TIE_BREAKING_INVARIANT_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3126"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3126_accepted_Pi_point_sources"], 0)
        self.assertEqual(cert["finite_Z12_orbit_obstruction_rows"], 10)
        self.assertEqual(cert["finite_orbit_selector_free_unique_points"], 0)
        self.assertEqual(cert["candidate_Omega_tie_sources"], 20)
        self.assertEqual(cert["invariant_law_rows"], 380)
        self.assertEqual(cert["orbit_witness_rows"], 280)
        self.assertEqual(cert["Pi_Lambda_Phi_Delta_Iota_Kappa_Tau_Xi_R_coupling_rows"], 240)
        self.assertEqual(cert["candidate_gate_rows"], 500)
        self.assertGreaterEqual(cert["promising_internal_tie_breakers"], 10)
        self.assertEqual(cert["accepted_Omega_tie_sources"], 0)

    def test_orbit_obstruction_and_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        orbit_rows = objs["finite_Z12_orbit_obstruction_rows"]
        self.assertTrue(all(not row["selector_free_unique_point_available"] for row in orbit_rows))
        candidates = objs["candidate_Omega_tie_sources"]
        self.assertTrue(any(row["candidate"] == "centered_phase_information_skew" and row["translation_safe"] and not row["inversion_safe"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "translation_averaged_phase_skew" and row["translation_safe"] and not row["nonzero_value_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "a_phi_weighted_orbit_moment" and row["C_phi_A_phi_preserved"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "observed_light_tie_breaker" and not row["strict_nadsoliton_data_only"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_oriented_tie_breaker" and not row["non_selector_premise"] for row in candidates))
        self.assertTrue(all(not row["accepted_Omega_tie_source"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["accepted_coupling_chain"] for row in objs["Pi_Lambda_Phi_Delta_Iota_Kappa_Tau_Xi_R_coupling_rows"]))

    def test_recommendation_and_docs(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["finite_Z12_orbit_obstruction_computed"])
        self.assertTrue(decision["positive_scoped_flags"]["internal_phase_information_tie_breakers_remain_promising_but_scoped"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("finite Z12 orbit calculation", decision["bounded_result"])
        self.assertIn("strict pointed-orbit source law Sigma_point", decision["next_honest_step"])
        self.assertIn("rather than a post-hoc tie breaker", decision["next_honest_step"])
        self.assertIn("P3127/S2077", MD.read_text(encoding="utf-8"))
        self.assertIn("P3127/S2077", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3127/S2077", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3127", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
