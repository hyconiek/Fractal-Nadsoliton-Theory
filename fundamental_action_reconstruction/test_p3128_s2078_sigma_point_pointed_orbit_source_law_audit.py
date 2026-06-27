import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3128_s2078_sigma_point_pointed_orbit_source_law_audit.py"
OUT = ROOT / "generated" / "p3128_s2078_sigma_point_pointed_orbit_source_law_audit.json"
MD = ROOT / "generated" / "p3128_s2078_sigma_point_pointed_orbit_source_law_audit.md"


class P3128SigmaPointPointedOrbitSourceLawAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3128_SIGMA_POINT_POINTED_ORBIT_SOURCE_LAW_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3127"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3127_accepted_Omega_tie_sources"], 0)
        self.assertEqual(cert["finite_signed_orbit_obstruction_rows"], 12)
        self.assertEqual(cert["finite_orbit_selector_free_signed_representatives"], 0)
        self.assertEqual(cert["candidate_Sigma_point_sources"], 21)
        self.assertEqual(cert["source_law_rows"], 441)
        self.assertEqual(cert["signed_orbit_witness_rows"], 336)
        self.assertEqual(cert["Omega_Pi_Lambda_Phi_Delta_Iota_Kappa_Tau_Xi_R_coupling_rows"], 273)
        self.assertEqual(cert["candidate_gate_rows"], 567)
        self.assertGreaterEqual(cert["promising_internal_signed_sources"], 10)
        self.assertEqual(cert["accepted_Sigma_point_sources"], 0)

    def test_signed_orbit_obstruction_and_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        orbit_rows = objs["finite_signed_orbit_obstruction_rows"]
        self.assertTrue(all(not row["selector_free_signed_representative_available"] for row in orbit_rows))
        self.assertTrue(all(set(row["inversion_pairs_signs"]) == {-1, 1} for row in orbit_rows))
        candidates = objs["candidate_Sigma_point_sources"]
        self.assertTrue(any(row["candidate"] == "translation_covariant_signed_density" and row["translation_covariant"] and not row["inversion_safe"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "aut_equivariant_signed_orbit_section" and row["aut_equivariant"] and not row["nonzero_signed_representative_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "a_phi_signed_cell_source" and row["C_phi_A_phi_preserved"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "omega_tie_repackaged_source" and not row["direct_source_not_posthoc_tie_breaker"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_signed_orbit_source" and not row["selector_free"] for row in candidates))
        self.assertTrue(all(not row["accepted_Sigma_point_source"] for row in objs["candidate_aggregate_certificate"]))
        self.assertTrue(all(not row["accepted_coupling_chain"] for row in objs["Omega_Pi_Lambda_Phi_Delta_Iota_Kappa_Tau_Xi_R_coupling_rows"]))

    def test_recommendation_and_docs(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["finite_signed_Z12_orbit_obstruction_computed"])
        self.assertTrue(decision["positive_scoped_flags"]["internal_phase_information_signed_sources_remain_promising_but_scoped"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("finite signed Z12 orbit calculation", decision["bounded_result"])
        self.assertIn("strict sign-and-origin generator object Gamma_SO", decision["next_honest_step"])
        self.assertIn("both a nonzero sign and a source-origin representative", decision["next_honest_step"])
        self.assertIn("P3128/S2078", MD.read_text(encoding="utf-8"))
        self.assertIn("P3128/S2078", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3128/S2078", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3128", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
