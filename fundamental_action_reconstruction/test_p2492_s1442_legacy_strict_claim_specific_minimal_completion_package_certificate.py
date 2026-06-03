from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2492_s1442_legacy_strict_claim_specific_minimal_completion_package_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2492LegacyStrictClaimSpecificMinimalCompletionPackageCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["legacy_strict_claim_specific_minimal_completion_package_certificate"]["theorem_export"]
        cls.cert = cls.theorem["claim_specific_minimal_completion_package_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2492")
        self.assertEqual(self.payload["stage_id"], "S1442")
        self.assertIn("MINIMAL_COMPLETION_PACKAGE", self.payload["status"])
        self.assertIn("NO_ATOM_EXPORT", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_package_counts(self) -> None:
        self.assertEqual(self.theorem["physical_claim_count"], 4)
        self.assertEqual(self.theorem["shared_core_atom_count"], 10)
        self.assertEqual(self.theorem["minimal_physical_claim_completion_size"], 11)
        self.assertEqual(self.theorem["maximal_physical_claim_completion_size"], 12)
        self.assertEqual(self.theorem["all_roles_minimal_completion_size"], 14)
        self.assertTrue(self.theorem["all_roles_package_equals_union_of_physical_claim_packages"])
        self.assertEqual(self.cert["bridge_shared_core_atom_count"], 8)
        self.assertEqual(self.cert["role_shared_core_atom_count"], 2)

    def test_claim_specific_rows(self) -> None:
        rows = {row["claim"]: row for row in self.cert["physical_claim_minimal_completion_rows"]}
        self.assertEqual(rows["legacy_weak_mixing_angle"]["minimal_current_completion_size"], 11)
        self.assertEqual(rows["legacy_inverse_alpha_em"]["minimal_current_completion_size"], 12)
        self.assertEqual(rows["legacy_beta_power_gravity_hierarchy"]["minimal_current_completion_size"], 12)
        self.assertEqual(rows["legacy_torsion_to_chi11_orientation"]["minimal_current_completion_size"], 12)
        self.assertEqual(rows["legacy_weak_mixing_angle"]["stage_counts"], {"bridge": 8, "role": 3})
        self.assertEqual(rows["legacy_inverse_alpha_em"]["stage_counts"], {"bridge": 8, "role": 4})
        self.assertIn("role::alpha_geo_strict_role_successor_theorem", rows["legacy_weak_mixing_angle"]["claim_specific_atoms"])
        self.assertIn("role::strict_nonlinear_hierarchy_successor_theorem", rows["legacy_beta_power_gravity_hierarchy"]["claim_specific_atoms"])
        self.assertIn("role::chi11_orientation_role_successor_theorem", rows["legacy_torsion_to_chi11_orientation"]["claim_specific_atoms"])

    def test_overlap_and_dependency_structure(self) -> None:
        self.assertTrue(self.theorem["weak_mixing_package_is_proper_subset_of_inverse_alpha_package"])
        self.assertTrue(self.theorem["role_transfer_and_ltotal_in_every_physical_package"])
        self.assertTrue(self.theorem["every_physical_package_contains_all_bridge_atoms"])
        self.assertEqual(set(self.cert["role_shared_core_atoms"]), {"role::role_transfer_audit_license", "role::role_bearing_ltotal_export"})
        self.assertEqual(set(self.cert["alpha_geo_successor_claims"]), {"legacy_inverse_alpha_em", "legacy_weak_mixing_angle"})
        self.assertEqual(set(self.cert["beta_tors_successor_claims"]), {"legacy_beta_power_gravity_hierarchy", "legacy_inverse_alpha_em", "legacy_torsion_to_chi11_orientation"})
        overlap = {tuple(row["claim_pair"]): row for row in self.cert["package_overlap_rows"]}
        weak_inverse = overlap[("legacy_inverse_alpha_em", "legacy_weak_mixing_angle")]
        self.assertTrue(weak_inverse["right_package_is_subset_of_left"])
        self.assertEqual(weak_inverse["jaccard_fraction"], "11/12")

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
        self.assertIn("P2492/S1442", MD.read_text(encoding="utf-8"))
        self.assertIn("P2492/S1442", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2492/S1442", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
