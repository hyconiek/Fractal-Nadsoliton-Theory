from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2509_s1459_strict_damping_rg_minimum_roughness_variational_wellposedness_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2509StrictDampingRGMinimumRoughnessVariationalWellposednessCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_rg_minimum_roughness_variational_wellposedness_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_rg_minimum_roughness_variational_wellposedness_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2509")
        self.assertEqual(self.payload["stage_id"], "S1459")
        self.assertIn("MINIMUM_ROUGHNESS_VARIATIONAL_WELLPOSEDNESS", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_THEOREM", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_variational_wellposedness(self) -> None:
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_damping_beta_eta_source")
        self.assertTrue(self.theorem["p2508_sobolev_node_coercivity_inherited"])
        self.assertTrue(self.theorem["candidate_matches_all_nodes"])
        self.assertEqual(float(self.theorem["max_abs_node_residual"]), 0.0)
        self.assertTrue(self.theorem["minimum_roughness_problem_wellposed_for_postulated_functional"])
        self.assertTrue(self.theorem["unique_minimizer_is_constant_flow_y0_delta_ell"])
        self.assertTrue(self.theorem["sample_tangent_all_node_vanishing"])
        self.assertTrue(self.theorem["sample_tangent_all_positive_energy"])

    def test_conditional_status_and_negative_controls(self) -> None:
        self.assertTrue(self.theorem["coercivity_is_selector_support_not_source_derivation"])
        self.assertTrue(self.theorem["roughness_action_still_postulated_not_derived"])
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
        self.assertIn("P2509/S1459", MD.read_text(encoding="utf-8"))
        self.assertIn("P2509/S1459", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2509/S1459", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
