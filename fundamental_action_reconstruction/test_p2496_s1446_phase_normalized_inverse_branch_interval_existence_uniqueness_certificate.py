from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2496_s1446_phase_normalized_inverse_branch_interval_existence_uniqueness_certificate import (
    MD,
    OUT,
    build_payload,
    write_markdown,
)


class P2496PhaseNormalizedInverseBranchIntervalExistenceUniquenessCertificateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        cls.theorem = cls.payload["phase_normalized_inverse_branch_interval_existence_uniqueness_certificate"]["theorem_export"]
        cls.cert = cls.theorem["inverse_branch_interval_existence_uniqueness_certificate"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2496")
        self.assertEqual(self.payload["stage_id"], "S1446")
        self.assertIn("INVERSE_BRANCH_INTERVAL_EXISTENCE_UNIQUENESS", self.payload["status"])
        self.assertIn("NO_GLOBAL_BRANCH_THEOREM", self.payload["status"])
        self.assertIn("NO_SOURCE_EXPORT", self.payload["status"])
        self.assertIn("NO_BRIDGE_ATOM", self.payload["status"])
        self.assertIn("NO_QW2191", self.payload["status"])
        self.assertIn("NO_ROLE_TRANSFER", self.payload["status"])
        self.assertIn("NO_TOE", self.payload["status"])

    def test_legacy_branch_derivative_interval(self) -> None:
        bounds = self.cert["legacy_derivative_interval_on_0_2"]
        self.assertEqual(self.cert["legacy_derivative_interval_backend"], "mpmath.iv")
        self.assertEqual(self.cert["legacy_derivative_interval_dps"], 100)
        self.assertLess(float(bounds[0]), 0.0)
        self.assertLess(float(bounds[1]), 0.0)
        self.assertTrue(self.theorem["legacy_derivative_interval_strictly_negative"])
        self.assertEqual(self.cert["legacy_branch_interval"], ["0.0", "2.0"])

    def test_sample_and_z12_targets_are_bracketed(self) -> None:
        self.assertTrue(self.theorem["all_sample_targets_inside_branch_range"])
        self.assertTrue(self.theorem["all_z12_targets_inside_branch_range"])
        self.assertEqual(len(self.cert["sample_branch_bracket_rows"]), 10)
        self.assertEqual(len(self.cert["z12_branch_bracket_rows"]), 12)
        self.assertTrue(all(row["inside_legacy_branch_range"] for row in self.cert["sample_branch_bracket_rows"]))
        self.assertTrue(all(row["inside_legacy_branch_range"] for row in self.cert["z12_branch_bracket_rows"]))
        z12_zero = self.cert["z12_branch_bracket_rows"][0]
        self.assertEqual(z12_zero["d"], "0")
        self.assertTrue(z12_zero["endpoint_case_at_left"])

    def test_margins_and_inherited_interval_enclosure(self) -> None:
        self.assertGreater(float(self.theorem["positive_right_margin_min"]), 0.55)
        self.assertGreater(float(self.theorem["nonzero_left_margin_min_excluding_d0"]), 3e-6)
        self.assertTrue(self.theorem["p2495_interval_enclosure_inherited"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))

    def test_negative_controls(self) -> None:
        self.assertFalse(self.theorem["formal_directed_rounding_backend_exported"])
        self.assertFalse(self.theorem["global_branch_theorem_exported"])
        self.assertFalse(self.theorem["curvature_dynamic_source_exported"])
        self.assertFalse(self.theorem["legacy_to_strict_bridge_atom_exported"])
        self.assertFalse(self.theorem["strict_compression_dynamic_source_exported"])
        self.assertFalse(self.theorem["selector_source_theorem_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_certificate"])
        self.assertFalse(self.theorem["role_transfer_licensed_by_this_certificate"])
        self.assertFalse(self.theorem["toe_closure_exported"])

    def test_rg_audit_and_documentation_markers(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        self.assertGreaterEqual(audit["patterns"]["new_packet"]["count"], 0)
        self.assertIn("P2496/S1446", MD.read_text(encoding="utf-8"))
        self.assertIn("P2496/S1446", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2496/S1446", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
