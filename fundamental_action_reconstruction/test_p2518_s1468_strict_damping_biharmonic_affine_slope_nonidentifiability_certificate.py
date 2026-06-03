from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2518_s1468_strict_damping_biharmonic_affine_slope_nonidentifiability_certificate import (
    MD,
    OUT,
    SLOPE_CANDIDATES,
    build_payload,
    write_markdown,
)


class P2518StrictDampingBiharmonicAffineSlopeNonidentifiabilityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_biharmonic_affine_slope_nonidentifiability_certificate"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_biharmonic_affine_slope_nonidentifiability_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2518")
        self.assertEqual(self.payload["stage_id"], "S1468")
        self.assertIn("BIHARMONIC_AFFINE_SLOPE_NONIDENTIFIABILITY", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_affine_slope_operator_nonidentifiability(self) -> None:
        self.assertTrue(self.theorem["p2517_axiom_boundary_inherited"])
        self.assertEqual(self.theorem["strict_running_slope_delta"], "4/5")
        self.assertEqual(self.theorem["affine_family_dimension"], 2)
        self.assertEqual(self.theorem["distinct_zero_energy_slope_count"], len(SLOPE_CANDIDATES))
        self.assertTrue(self.theorem["all_audited_slopes_operator_indistinguishable"])
        self.assertTrue(self.theorem["strict_delta_is_only_one_member_of_zero_energy_family"])
        self.assertTrue(all(row["m2_operator_signature_observables_all_zero"] for row in self.cert["affine_operator_observable_rows"]))

    def test_finite_rank_audit(self) -> None:
        audit = self.cert["finite_operator_rank_audit"]
        self.assertTrue(audit["rank_identities_pass"])
        self.assertTrue(audit["D2_kernel_dimension_is_two"])
        self.assertTrue(audit["D4_kernel_dimension_contains_affine_but_is_larger"])
        row_by_operator = {row["operator"]: row for row in audit["rows"]}
        self.assertEqual(row_by_operator["D^2"]["kernel_dimension"], 2)
        self.assertEqual(row_by_operator["D^4"]["kernel_dimension"], 4)

    def test_negative_controls(self) -> None:
        self.assertTrue(self.theorem["operator_signature_numeric_key_separation_exported"])
        self.assertFalse(self.theorem["m2_operator_signature_source_exported"])
        self.assertFalse(self.theorem["beta_eta_numeric_source_exported"])
        self.assertFalse(self.theorem["strict_damping_beta_eta_source_exported"])
        self.assertFalse(self.theorem["damping_compression_bridge_component_ready"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["selector_closure_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2518/S1468", MD.read_text(encoding="utf-8"))
        self.assertIn("P2518/S1468", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2518/S1468", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
