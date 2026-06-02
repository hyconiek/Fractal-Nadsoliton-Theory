import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2443_s1393_strict_moment_supplemental_constraint_rank_lift_certificate.py"
OUT = ROOT / "generated" / "p2443_s1393_strict_moment_supplemental_constraint_rank_lift_certificate.json"
MD = ROOT / "generated" / "p2443_s1393_strict_moment_supplemental_constraint_rank_lift_certificate.md"
P2442 = ROOT / "generated" / "p2442_s1392_strict_moment_coefficient_local_identifiability_nullspace_certificate.json"


class P2443StrictMomentSupplementalConstraintRankLiftCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2442.exists():
            subprocess.run([sys.executable, str(ROOT / "p2442_s1392_strict_moment_coefficient_local_identifiability_nullspace_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["strict_moment_supplemental_constraint_rank_lift_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_outputs_and_candidate_count(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2443")
        self.assertEqual(self.payload["stage_id"], "S1393")
        self.assertEqual(self.payload["status"], "PASS_STRICT_MOMENT_SUPPLEMENTAL_CONSTRAINT_RANK_LIFT_FRONTIER_NO_SOURCE_THEOREM")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["inherited_base_rank"], 3)
        self.assertEqual(self.theorem["inherited_nullspace_dimension"], 1)
        self.assertEqual(self.theorem["candidate_count"], 12)

    def test_rank_lift_frontier(self) -> None:
        self.assertEqual(self.theorem["rank_lifting_candidate_count"], 12)
        self.assertEqual(self.theorem["minimal_rank_lift_antichain_size"], 1)
        self.assertEqual(len(self.theorem["minimal_rank_lift_singletons"]), 12)
        self.assertTrue(self.theorem["all_raw_moment_candidates_rank_lift"])
        self.assertTrue(self.theorem["all_kernel_sample_candidates_rank_lift"])
        self.assertTrue(all(row["augmented_rank"] == 4 for row in self.theorem["candidate_rows"]))

    def test_hard_limits(self) -> None:
        self.assertFalse(self.theorem["strict_observable_source_constraint_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_kernel_to_coefficient_map_theorem_exported"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_docs(self) -> None:
        self.assertIn("supplemental-constraint rank-lift certificate", MD.read_text(encoding="utf-8"))
        self.assertIn("P2443/S1393", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2443/S1393", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
