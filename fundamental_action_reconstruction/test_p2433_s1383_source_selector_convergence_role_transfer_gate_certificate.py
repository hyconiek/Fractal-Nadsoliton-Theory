#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2433_s1383_source_selector_convergence_role_transfer_gate_certificate.py"
OUT = ROOT / "generated" / "p2433_s1383_source_selector_convergence_role_transfer_gate_certificate.json"
MD = ROOT / "generated" / "p2433_s1383_source_selector_convergence_role_transfer_gate_certificate.md"
P2432 = ROOT / "generated" / "p2432_s1382_post_antichain_branch_residual_transition_certificate.json"


class P2433SourceSelectorConvergenceRoleTransferGateCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not P2432.exists():
            subprocess.run([sys.executable, str(ROOT / "p2432_s1382_post_antichain_branch_residual_transition_certificate.py")], check=True)
        subprocess.run([sys.executable, str(SCRIPT)], check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cert = cls.payload["source_selector_convergence_role_transfer_gate_certificate"]
        cls.theorem = cls.cert["theorem_export"]

    def test_identity_and_outputs(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2433")
        self.assertEqual(self.payload["stage_id"], "S1383")
        self.assertEqual(
            self.payload["status"],
            "PASS_SOURCE_SELECTOR_CONVERGENCE_ROLE_TRANSFER_NEXT_NO_GATE_DISCHARGE_NO_CLOSURE",
        )
        self.assertTrue(OUT.exists())
        self.assertTrue(MD.exists())

    def test_convergence_state(self) -> None:
        self.assertEqual(self.theorem["convergence_state_count"], 2)
        self.assertTrue(self.theorem["convergence_states_have_same_gate_set"])
        self.assertEqual(
            self.theorem["readiness_after_convergence"],
            {
                "bridge_source_ready": True,
                "selector_source_ready": True,
                "role_transfer_ready": False,
                "role_bearing_ltotal_ready": False,
                "toe_ready": False,
            },
        )
        self.assertEqual(
            self.theorem["remaining_gates_after_convergence"],
            ["role_transfer_audit_license", "role_bearing_ltotal_export"],
        )

    def test_role_transfer_next_not_closure(self) -> None:
        self.assertTrue(self.theorem["role_transfer_admissible_after_source_selector_convergence"])
        self.assertFalse(self.theorem["role_bearing_ltotal_admissible_after_source_selector_convergence"])
        self.assertFalse(self.theorem["toe_ready_after_source_selector_convergence"])
        self.assertTrue(self.theorem["p2432_source_first_next_inherited"])
        self.assertTrue(self.theorem["p2432_selector_pair_next_inherited"])
        self.assertTrue(self.theorem["p2432_role_transfer_not_after_first_branches_inherited"])

    def test_hard_limits_and_docs(self) -> None:
        self.assertFalse(self.theorem["source_obligation_discharge_exported_by_this_certificate"])
        self.assertFalse(self.theorem["chi11_source_exported_by_this_certificate"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_transfer_licensed"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertIn("source-selector convergence", MD.read_text(encoding="utf-8"))
        self.assertIn("P2433/S1383", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2433/S1383", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
