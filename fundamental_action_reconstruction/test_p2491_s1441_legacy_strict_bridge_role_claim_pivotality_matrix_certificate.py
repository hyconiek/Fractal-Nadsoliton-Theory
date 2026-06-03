from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2491_s1441_legacy_strict_bridge_role_claim_pivotality_matrix_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2491LegacyStrictBridgeRoleClaimPivotalityMatrixCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["legacy_strict_bridge_role_claim_pivotality_matrix_certificate"]["theorem_export"]
        cls.matrix = cls.theorem["bridge_role_claim_pivotality_matrix"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2491")
        self.assertEqual(self.payload["stage_id"], "S1441")
        self.assertIn("PIVOTALITY_MATRIX", self.payload["status"])
        self.assertIn("NO_ATOM_EXPORT", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_matrix_counts_and_current_state(self) -> None:
        self.assertEqual(self.theorem["combined_atom_count"], 14)
        self.assertEqual(self.theorem["assignment_count"], 16384)
        self.assertEqual(self.theorem["one_atom_flip_case_count"], 114688)
        self.assertEqual(self.theorem["claim_count_including_all_roles_target"], 5)
        self.assertEqual(self.theorem["current_ready_claim_count"], 0)
        self.assertEqual(len(self.matrix["matrix_rows"]), 14)
        self.assertEqual(len(self.matrix["claim_requirements"]), 5)

    def test_pivotality_structure(self) -> None:
        self.assertTrue(self.theorem["all_bridge_atoms_have_equal_physical_role_pivotal_total"])
        self.assertTrue(self.theorem["role_transfer_and_ltotal_are_top_role_stage_pivotal_atoms"])
        self.assertTrue(self.theorem["every_atom_is_pivotal_for_all_roles_target_once"])
        self.assertEqual(self.theorem["top_physical_role_claim_pivotal_total"], 20)
        self.assertEqual(set(self.theorem["top_role_stage_physical_role_claim_pivotal_atoms"]), {"role::role_transfer_audit_license", "role::role_bearing_ltotal_export"})
        bridge_rows = [row for row in self.matrix["matrix_rows"] if row["stage"] == "bridge"]
        self.assertEqual({row["physical_role_claim_pivotal_total"] for row in bridge_rows}, {20})
        role_rows = {row["atom"]: row for row in self.matrix["matrix_rows"] if row["stage"] == "role"}
        self.assertEqual(role_rows["role::role_transfer_audit_license"]["physical_role_claim_pivotal_total"], 20)
        self.assertEqual(role_rows["role::role_bearing_ltotal_export"]["physical_role_claim_pivotal_total"], 20)
        self.assertEqual(role_rows["role::beta_tors_strict_role_successor_theorem"]["physical_role_claim_pivotal_total"], 12)

    def test_claim_specific_counts(self) -> None:
        role_rows = {row["atom"]: row for row in self.matrix["matrix_rows"]}
        alpha = role_rows["role::alpha_geo_strict_role_successor_theorem"]["claim_pivotal_counts"]
        self.assertEqual(alpha["legacy_weak_mixing_angle"]["pivotal_count"], 8)
        self.assertEqual(alpha["legacy_inverse_alpha_em"]["pivotal_count"], 4)
        self.assertEqual(alpha["legacy_beta_power_gravity_hierarchy"]["pivotal_count"], 0)
        self.assertEqual(alpha["legacy_torsion_to_chi11_orientation"]["pivotal_count"], 0)
        chi11 = role_rows["role::chi11_orientation_role_successor_theorem"]["claim_pivotal_counts"]
        self.assertEqual(chi11["legacy_torsion_to_chi11_orientation"]["pivotal_count"], 4)
        self.assertEqual(chi11["all_four_audited_legacy_role_claims"]["pivotal_count"], 1)

    def test_gatekeepers_and_negative_controls(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertFalse(self.theorem["source_obligation_discharge_exported_by_this_certificate"])
        self.assertFalse(self.theorem["selector_source_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_transfer_licensed_by_this_certificate"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["legacy_role_transfer_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2491/S1441", MD.read_text(encoding="utf-8"))
        self.assertIn("P2491/S1441", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2491/S1441", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
