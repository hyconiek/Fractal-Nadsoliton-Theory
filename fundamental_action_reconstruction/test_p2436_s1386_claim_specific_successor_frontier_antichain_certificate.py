import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2436_s1386_claim_specific_successor_frontier_antichain_certificate.py"
OUT = ROOT / "generated" / "p2436_s1386_claim_specific_successor_frontier_antichain_certificate.json"
MD = ROOT / "generated" / "p2436_s1386_claim_specific_successor_frontier_antichain_certificate.md"
P2435 = ROOT / "generated" / "p2435_s1385_legacy_role_claim_implication_separability_certificate.json"


class P2436ClaimSpecificSuccessorFrontierAntichainCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2435.exists():
            subprocess.run([sys.executable, str(ROOT / "p2435_s1385_legacy_role_claim_implication_separability_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["claim_specific_successor_frontier_antichain_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_outputs_and_counts(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2436")
        self.assertEqual(self.payload["stage_id"], "S1386")
        self.assertEqual(
            self.payload["status"],
            "PASS_CLAIM_SPECIFIC_SUCCESSOR_FRONTIER_ANTICHAIN_NO_ROLE_TRANSFER_NO_CLOSURE",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())
        self.assertEqual(self.theorem["claim_specific_successor_count"], 4)
        self.assertEqual(self.theorem["successor_lattice_assignment_count"], 16)

    def test_frontier_distribution_and_singletons(self) -> None:
        self.assertEqual(self.theorem["current_successor_mask"], 0)
        self.assertEqual(self.theorem["current_ready_role_claim_count"], 0)
        self.assertEqual(self.theorem["ready_role_count_distribution"], {"0": 5, "1": 6, "2": 2, "3": 2, "4": 1})
        self.assertEqual(
            self.theorem["singleton_unlocking_successors_from_empty"],
            ["alpha_geo_strict_role_successor_theorem"],
        )
        self.assertEqual(len(self.theorem["singleton_nonunlocking_successors_from_empty"]), 3)

    def test_minimal_antichain_and_inheritance(self) -> None:
        self.assertEqual(self.theorem["minimal_successor_size_by_role_claim"]["legacy_weak_mixing_angle"], 1)
        for role_id in [
            "legacy_inverse_alpha_em",
            "legacy_beta_power_gravity_hierarchy",
            "legacy_torsion_to_chi11_orientation",
        ]:
            self.assertEqual(self.theorem["minimal_successor_size_by_role_claim"][role_id], 2)
        self.assertEqual(self.theorem["minimal_all_role_successor_antichain_size"], 4)
        self.assertEqual(self.theorem["all_role_ready_masks"], [15])
        self.assertTrue(self.theorem["p2435_rank_four_inherited"])
        self.assertTrue(self.theorem["p2435_single_nontrivial_implication_inherited"])

    def test_hard_limits_and_docs(self) -> None:
        self.assertTrue(self.theorem["common_gates_are_assumptions_not_exports"])
        self.assertFalse(self.theorem["legacy_role_claim_transferred_by_this_certificate"])
        self.assertFalse(self.theorem["common_role_transfer_gates_exported_by_this_certificate"])
        self.assertFalse(self.theorem["claim_specific_successor_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("claim-specific successor frontier", MD.read_text(encoding="utf-8"))
        self.assertIn("P2436/S1386", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2436/S1386", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
