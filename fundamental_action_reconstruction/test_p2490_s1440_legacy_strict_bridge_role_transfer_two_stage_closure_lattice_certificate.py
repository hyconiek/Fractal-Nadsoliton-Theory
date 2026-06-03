from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2490_s1440_legacy_strict_bridge_role_transfer_two_stage_closure_lattice_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2490LegacyStrictBridgeRoleTransferTwoStageClosureLatticeCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["legacy_strict_bridge_role_transfer_two_stage_closure_lattice_certificate"]["theorem_export"]
        cls.lattice = cls.theorem["two_stage_closure_lattice"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2490")
        self.assertEqual(self.payload["stage_id"], "S1440")
        self.assertIn("TWO_STAGE", self.payload["status"])
        self.assertIn("NO_BRIDGE_ATOM_EXPORT", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_combined_lattice_counts_and_current_state(self) -> None:
        self.assertEqual(self.theorem["bridge_atom_count"], 8)
        self.assertEqual(self.theorem["role_obligation_count"], 6)
        self.assertEqual(self.theorem["combined_atom_count"], 14)
        self.assertEqual(self.theorem["combined_assignment_count"], 16384)
        self.assertEqual(self.theorem["end_to_end_ready_assignment_count"], 1)
        self.assertEqual(self.theorem["end_to_end_ready_fraction"], "1/16384")
        self.assertFalse(self.theorem["current_bridge_ready"])
        self.assertEqual(self.theorem["current_ready_role_claim_count"], 0)
        self.assertFalse(self.theorem["current_end_to_end_ready"])
        self.assertEqual(self.theorem["minimum_new_atoms_needed_from_current_for_end_to_end_all_role_closure"], 14)

    def test_nearest_misses_and_distributions(self) -> None:
        self.assertEqual(self.theorem["nearest_miss_count"], 14)
        self.assertEqual(len(self.lattice["nearest_miss_rows"]), 14)
        for row in self.lattice["nearest_miss_rows"]:
            self.assertEqual(len(row["missing_bridge_atoms"]) + len(row["missing_role_obligations"]), 1)
        self.assertEqual(self.lattice["bridge_ready_assignment_count_in_combined_lattice"], 64)
        self.assertEqual(self.lattice["bridge_blocked_assignment_count_in_combined_lattice"], 16320)
        self.assertEqual(self.lattice["role_count_distribution_when_bridge_ready"], {"0": 53, "1": 6, "2": 2, "3": 2, "4": 1})
        self.assertEqual(self.lattice["role_count_distribution_when_bridge_blocked"], {"0": 13515, "1": 1530, "2": 510, "3": 510, "4": 255})

    def test_separate_gate_and_frontier_loads(self) -> None:
        self.assertTrue(self.theorem["bridge_completion_and_role_transfer_separate_gates"])
        self.assertTrue(self.theorem["completed_bridge_alone_would_not_transfer_legacy_roles"])
        top_atoms = [row["atom"] for row in self.theorem["top_frontier_atom_load_rows"]]
        self.assertIn("role_transfer_audit_license", top_atoms)
        self.assertIn("role_bearing_ltotal_export", top_atoms)
        self.assertIn("chi11_selector_source_theorem", top_atoms)
        self.assertEqual(self.lattice["p2411_bridge_ready_true_masks_inherited"], [255])
        self.assertEqual(self.lattice["p2434_all_roles_ready_masks_inherited"], [63])

    def test_gatekeepers_and_negative_controls(self) -> None:
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.theorem["current_legacy_role_claims_transferred_by_this_certificate"], 0)
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
        self.assertIn("P2490/S1440", MD.read_text(encoding="utf-8"))
        self.assertIn("P2490/S1440", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2490/S1440", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
