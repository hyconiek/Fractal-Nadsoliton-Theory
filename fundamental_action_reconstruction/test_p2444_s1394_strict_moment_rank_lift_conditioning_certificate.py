import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2444_s1394_strict_moment_rank_lift_conditioning_certificate.py"
OUT = ROOT / "generated" / "p2444_s1394_strict_moment_rank_lift_conditioning_certificate.json"
MD = ROOT / "generated" / "p2444_s1394_strict_moment_rank_lift_conditioning_certificate.md"
P2443 = ROOT / "generated" / "p2443_s1393_strict_moment_supplemental_constraint_rank_lift_certificate.json"


class P2444StrictMomentRankLiftConditioningCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2443.exists():
            subprocess.run([sys.executable, str(ROOT / "p2443_s1393_strict_moment_supplemental_constraint_rank_lift_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_moment_rank_lift_conditioning_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_outputs_and_conditioning_size(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2444")
        self.assertEqual(self.payload["stage_id"], "S1394")
        self.assertEqual(self.payload["status"], "PASS_STRICT_MOMENT_RANK_LIFT_CONDITIONING_FRONTIER_NO_SOURCE_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["inherited_candidate_count"], 12)
        self.assertEqual(len(self.theorem["conditioned_candidate_rows_descending"]), 12)

    def test_conditioning_frontier(self) -> None:
        self.assertEqual(self.theorem["best_conditioned_candidate_id"], "K_at_d_1")
        self.assertGreater(self.theorem["best_conditioned_normalized_volume"], self.theorem["conditioning_threshold_normalized_volume"])
        self.assertEqual(self.theorem["robust_rank_lift_candidate_count"], 6)
        self.assertIn("M0", self.theorem["robust_rank_lift_candidate_ids"])
        self.assertEqual(self.theorem["weakest_conditioned_candidate_id"], "K_at_d_20")
        self.assertLess(self.theorem["weakest_conditioned_normalized_volume"], self.theorem["conditioning_threshold_normalized_volume"])

    def test_hard_limits(self) -> None:
        self.assertFalse(self.theorem["strict_observable_source_constraint_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("rank-lift conditioning certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2444/S1394", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2444/S1394", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
