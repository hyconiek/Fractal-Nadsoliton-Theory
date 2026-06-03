from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2512_s1462_strict_damping_rg_quadratic_source_admissibility_audit import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2512StrictDampingRGQuadraticSourceAdmissibilityAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["strict_damping_rg_quadratic_source_admissibility_audit"]["theorem_export"]
        cls.cert = cls.theorem["strict_damping_rg_quadratic_source_admissibility_audit"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2512")
        self.assertEqual(self.payload["stage_id"], "S1462")
        self.assertIn("QUADRATIC_SOURCE_ADMISSIBILITY", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_THEOREM", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_source_admissibility_obstruction_and_ambiguity(self) -> None:
        self.assertEqual(self.theorem["frontier_atom_under_attack"], "strict_damping_beta_eta_source")
        self.assertTrue(self.theorem["p2511_spline_collapse_inherited"])
        self.assertTrue(self.theorem["mass_term_obstruction_witness_exists"])
        self.assertTrue(self.theorem["derivative_only_quadratics_stationary_on_audited_tangent_basis"])
        self.assertTrue(self.theorem["derivative_only_rows_admitted"])
        self.assertTrue(self.theorem["mass_term_rows_rejected_without_extra_forcing"])
        self.assertTrue(self.theorem["mass_term_obstruction_theorem_exported_for_unforced_quadratic_sources"])
        self.assertTrue(self.theorem["derivative_only_source_ambiguity_identified"])
        self.assertTrue(self.theorem["roughness_order_not_uniquely_sourced_by_stationarity"])
        witness = self.cert["stationarity_witness_audit"]
        self.assertTrue(witness["all_basis_node_vanishing"])
        self.assertGreater(witness["nonzero_mass_obstruction_row_count"], 0)
        self.assertGreater(float(witness["max_abs_mass_variation_moment"]), 1e-12)
        self.assertLess(float(witness["max_abs_first_derivative_variation"]), 1e-80)
        self.assertLess(float(witness["max_abs_roughness_variation"]), 1e-80)

    def test_conditional_status_and_negative_controls(self) -> None:
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
        self.assertIn("P2512/S1462", MD.read_text(encoding="utf-8"))
        self.assertIn("P2512/S1462", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2512/S1462", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
