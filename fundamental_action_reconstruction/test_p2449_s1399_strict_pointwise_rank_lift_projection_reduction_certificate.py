import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.py"
OUT = ROOT / "generated" / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json"
MD = ROOT / "generated" / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.md"
P2448 = ROOT / "generated" / "p2448_s1398_strict_pointwise_rank_lift_global_stationary_census_certificate.json"


class P2449StrictPointwiseRankLiftProjectionReductionCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2448.exists():
            subprocess.run([sys.executable, str(ROOT / "p2448_s1398_strict_pointwise_rank_lift_global_stationary_census_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_pointwise_rank_lift_projection_reduction_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_outputs_and_projection_identity(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2449")
        self.assertEqual(self.payload["stage_id"], "S1399")
        self.assertEqual(self.payload["status"], "PASS_STRICT_POINTWISE_PROJECTION_REDUCTION_MATCHES_CENSUS_NO_SELECTOR_SOURCE_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertGreater(self.theorem["projection_vector_norm"], 1.0e-6)
        self.assertLess(self.theorem["maximum_projection_identity_error"], 1.0e-12)

    def test_stationary_factorization_matches_p2448(self) -> None:
        self.assertEqual(len(self.theorem["zero_projection_roots"]), 2)
        self.assertEqual(len(self.theorem["stationary_factor_roots"]), 1)
        self.assertEqual(self.theorem["reconstructed_root_count"], 3)
        self.assertTrue(self.theorem["all_reconstructed_roots_match_p2448"])
        self.assertTrue(self.theorem["best_reconstructed_matches_p2448_best"])
        best = self.theorem["best_reconstructed_root"]
        self.assertEqual(best["root_family"], "stationary_factor")
        self.assertGreater(best["root_d"], 0.785)
        self.assertLess(best["root_d"], 0.786)
        self.assertGreater(best["normalized_rank_lift_volume"], 0.022)

    def test_hard_limits(self) -> None:
        self.assertFalse(self.theorem["exact_interval_root_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["pointwise_coordinate_selector_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_observable_source_constraint_exported_by_this_certificate"])
        self.assertFalse(self.theorem["gauge_slice_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("projection-reduction certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2449/S1399", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2449/S1399", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
