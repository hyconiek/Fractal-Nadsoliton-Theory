from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2506StrictDampingRGMinimumRoughnessSelectorCandidateCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_rg_minimum_roughness_selector_candidate_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_rg_minimum_roughness_selector_candidate_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2506")
        self.assertEqual(self.payload["stage_id"], "S1456")
        self.assertIn("MINIMUM_ROUGHNESS_SELECTOR_CANDIDATE", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_THEOREM", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_selector_candidate_and_nullspace_energy(self) -> None:
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_damping_beta_eta_source")
        self.assertTrue(self.theorem["p2505_nullspace_inherited"])
        self.assertTrue(self.theorem["constant_candidate_energy_zero"])
        self.assertTrue(self.theorem["selector_candidate_unique_if_roughness_action_is_admitted"])
        self.assertGreater(float(self.theorem["nullspace_scaled_roughness_energy"]), 0.0)
        self.assertTrue(self.theorem["nullspace_scaled_energy_positive"])
        self.assertGreater(float(self.theorem["max_abs_epsilon_R_second_midpoint"]), 0.0)
        self.assertTrue(self.theorem["constant_flow_selected_if_selector_postulated"])

    def test_conditional_status_and_negative_controls(self) -> None:
        self.assertTrue(self.theorem["selector_is_postulated_not_derived"])
        self.assertTrue(self.theorem["finite_node_nullspace_blocker_addressed_conditionally"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["strict_dynamical_source_for_A_P_D_exported"])
        self.assertFalse(self.theorem["strict_phase_frequency_source_exported"])
        self.assertFalse(self.theorem["bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["selector_closure_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2506/S1456", MD.read_text(encoding="utf-8"))
        self.assertIn("P2506/S1456", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2506/S1456", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
