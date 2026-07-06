import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
SCRIPT = ROOT / "p3129_s2079_gamma_so_sign_origin_generator_audit.py"
OUT = ROOT / "generated" / "p3129_s2079_gamma_so_sign_origin_generator_audit.json"
MD = ROOT / "generated" / "p3129_s2079_gamma_so_sign_origin_generator_audit.md"


class P3129GammaSOSignOriginGeneratorAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True, stdout=subprocess.PIPE, text=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_status_and_predecessor_hash(self):
        self.assertEqual(self.payload["status"], "P3129_GAMMA_SO_SIGN_ORIGIN_GENERATOR_BOUNDED_NO_GO")
        self.assertTrue(self.payload["input_hashes"]["P3128"])

    def test_certificate_counts(self):
        cert = self.payload["finite_certificate"]
        self.assertEqual(cert["content_grep_lanes"], 4)
        self.assertGreater(cert["content_grep_hits"], 0)
        self.assertEqual(cert["p3128_accepted_Sigma_point_sources"], 0)
        self.assertEqual(cert["finite_joint_sign_origin_obstruction_rows"], 12)
        self.assertEqual(cert["selector_free_joint_sign_origin_rows"], 0)
        self.assertEqual(cert["candidate_Gamma_SO_generators"], 16)
        self.assertEqual(cert["generator_law_rows"], 256)
        self.assertEqual(cert["symmetry_witness_rows"], 208)
        self.assertEqual(cert["candidate_gate_rows"], 336)
        self.assertGreaterEqual(cert["promising_internal_Gamma_SO_candidates"], 8)
        self.assertEqual(cert["accepted_Gamma_SO_generators"], 0)
        self.assertGreater(cert["A_phi"], 0)

    def test_finite_obstruction_and_candidate_rejections(self):
        objs = self.payload["constructed_theoretical_objects"]
        finite_rows = objs["finite_joint_sign_origin_obstruction_rows"]
        self.assertTrue(all(not row["gamma_so_selector_free_available"] for row in finite_rows))
        self.assertTrue(all(set(row["sign_values_seen"]) == {-1, 1} for row in finite_rows))
        self.assertTrue(all(len(row["origin_values_seen"]) > 1 for row in finite_rows))
        candidates = objs["candidate_Gamma_SO_generators"]
        self.assertTrue(any(row["candidate"] == "chiral_bispectrum_argmax_origin" and row["nonzero_sign_exported"] and row["source_origin_representative_exported"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "category_initial_pointed_sign" and row["inversion_safe"] and not row["Sigma_point_retest_unlocked"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "selector_chosen_gamma_so" and not row["selector_free"] for row in candidates))
        self.assertTrue(any(row["candidate"] == "least_action_gamma_so" and not row["not_imported_dynamics"] for row in candidates))
        self.assertTrue(all(not row["accepted_Gamma_SO_generator"] for row in objs["candidate_aggregate_certificate"]))

    def test_recommendation_and_docs(self):
        decision = self.payload["decision"]
        self.assertTrue(decision["positive_scoped_flags"]["finite_joint_sign_origin_obstruction_computed"])
        self.assertTrue(decision["positive_scoped_flags"]["internal_sign_origin_candidates_remain_promising_but_scoped"])
        self.assertTrue(all(value is False for value in decision["negative_export_flags"].values()))
        self.assertIn("finite joint sign-origin calculation", decision["bounded_result"])
        self.assertIn("Theta_TO", decision["next_honest_step"])
        self.assertIn("translation-origin quotient", decision["next_honest_step"])
        self.assertIn("P3129/S2079", MD.read_text(encoding="utf-8"))
        self.assertIn("P3129/S2079", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P3129/S2079", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))
        self.assertIn("P3129", (REPO / "AGENTS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
