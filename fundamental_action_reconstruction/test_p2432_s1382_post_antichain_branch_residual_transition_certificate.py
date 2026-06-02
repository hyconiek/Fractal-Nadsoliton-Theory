#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2432_s1382_post_antichain_branch_residual_transition_certificate.py"
OUT = ROOT / "generated" / "p2432_s1382_post_antichain_branch_residual_transition_certificate.json"
MD = ROOT / "generated" / "p2432_s1382_post_antichain_branch_residual_transition_certificate.md"
P2431 = ROOT / "generated" / "p2431_s1381_admissible_next_theorem_target_antichain_certificate.json"


class P2432PostAntichainBranchResidualTransitionCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2431.exists():
            subprocess.run([sys.executable, str(ROOT / "p2431_s1381_admissible_next_theorem_target_antichain_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["post_antichain_branch_residual_transition_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2432")
        self.assertEqual(self.payload["stage_id"], "S1382")
        self.assertEqual(self.payload["status"], "PASS_POST_ANTICHAIN_BRANCH_TRANSITION_NO_GATE_DISCHARGE_NO_CLOSURE")
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_branch_transitions(self) -> None:
        self.assertEqual(self.theorem["branch_count"], 2)
        self.assertEqual(self.theorem["source_first_remaining_gate_count"], 4)
        self.assertEqual(self.theorem["selector_pair_first_remaining_gate_count"], 3)
        self.assertEqual(
            self.theorem["source_first_next_minimal_readiness_candidates"],
            [["chi11_source_export", "qw2191_selector_discharge"]],
        )
        self.assertEqual(
            self.theorem["selector_pair_first_next_minimal_readiness_candidates"],
            [["source_obligation_discharge"], ["source_obligation_discharge", "role_transfer_audit_license"]],
        )

    def test_role_blocks_and_inheritance(self) -> None:
        self.assertFalse(self.theorem["role_transfer_admissible_after_source_first"])
        self.assertFalse(self.theorem["role_transfer_admissible_after_selector_pair_first"])
        self.assertFalse(self.theorem["role_bearing_ltotal_admissible_after_either_first_branch"])
        self.assertTrue(self.theorem["p2431_antichain_inherited"])
        self.assertTrue(self.theorem["p2431_role_transfer_not_first_inherited"])
        self.assertTrue(self.theorem["p2431_ltotal_not_first_inherited"])
        self.assertTrue(self.theorem["p2430_global_minimal_cover_inherited"])

    def test_hard_limits_and_docs(self) -> None:
        self.assertFalse(self.theorem["source_obligation_discharge_exported"])
        self.assertFalse(self.theorem["chi11_source_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("post-antichain branch", MD.read_text(encoding="utf-8"))
        self.assertIn("P2432/S1382", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2432/S1382", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
