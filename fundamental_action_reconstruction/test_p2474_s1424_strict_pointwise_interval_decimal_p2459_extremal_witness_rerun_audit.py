#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import p2474_s1424_strict_pointwise_interval_decimal_p2459_extremal_witness_rerun_audit as p2474

OUT = ROOT / "generated" / "p2474_s1424_strict_pointwise_interval_decimal_p2459_extremal_witness_rerun_audit.json"
MD = ROOT / "generated" / "p2474_s1424_strict_pointwise_interval_decimal_p2459_extremal_witness_rerun_audit.md"


class P2474StrictPointwiseIntervalDecimalP2459ExtremalWitnessRerunAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not OUT.exists():
            p2474.main()
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))
        cls.theorem = cls.payload["strict_pointwise_interval_decimal_p2459_extremal_witness_rerun_audit"]["theorem_export"]

    def test_identity_and_status(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2474")
        self.assertEqual(self.payload["stage_id"], "S1424")
        self.assertIn("P2459_EXTREMAL_WITNESS_RERUN_AUDIT", self.payload["status"])
        self.assertIn("NO_DIRECTED_ROUNDING", self.payload["status"])

    def test_witness_rerun_results(self) -> None:
        self.assertEqual(self.theorem["total_fresh_decimal_taylor_witness_rerun_count"], 28)
        self.assertTrue(self.theorem["all_fresh_witness_groups_match_stored"])
        self.assertTrue(self.theorem["all_fresh_witness_groups_exclude_zero_with_positive_separation"])
        expected_counts = {"P2469/S1419": 6, "P2470/S1420": 6, "P2472/S1422": 16}
        for group in self.theorem["witness_groups"]:
            self.assertEqual(group["fresh_witness_count"], expected_counts[group["source_packet"]])
            self.assertTrue(group["all_fresh_witnesses_match_stored"])
            self.assertTrue(group["all_fresh_witnesses_exclude_zero_with_positive_separation"])
            for row in group["fresh_witness_rows"]:
                self.assertTrue(row["separation_matches_stored"])
                self.assertTrue(row["zero_exclusion_matches_stored"])
                self.assertTrue(row["fresh_decimal_separation_positive"])

    def test_chain_counts_and_inherited_binding(self) -> None:
        self.assertEqual(self.theorem["p2466_tail_count_inherited_from_p2471"], 6361)
        self.assertEqual(self.theorem["p2469_full_opposite_tail_count_inherited"], 45165)
        self.assertEqual(self.theorem["p2470_remaining_non_tail_count_inherited"], 48320)
        self.assertEqual(self.theorem["p2471_p2459_universe_count_inherited"], 99846)
        self.assertEqual(self.theorem["finite_chain_sum_p2466_p2469_p2470"], 99846)
        self.assertTrue(self.theorem["finite_chain_sum_matches_p2471_universe"])
        self.assertEqual(self.theorem["p2471_missing_cells_inherited"], 0)
        self.assertEqual(self.theorem["p2471_extra_cells_inherited"], 0)
        self.assertTrue(self.theorem["p2471_disjoint_complete_inherited"])
        self.assertTrue(self.theorem["p2473_fingerprint_binding_pass_inherited"])

    def test_gatekeepers_and_negative_controls(self) -> None:
        for name, value in self.payload["gatekeeper_checks"].items():
            self.assertTrue(value, name)
        self.assertTrue(self.theorem["finite_replay_chain_extremal_witness_rerun_audit_exported"])
        self.assertFalse(self.theorem["directed_rounding_interval_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["symbolic_root_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["analytic_monotonicity_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["global_continuum_root_exclusion_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["pointwise_coordinate_selector_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_observable_source_constraint_exported_by_this_certificate"])
        self.assertFalse(self.theorem["gauge_slice_theorem_exported_by_this_certificate"])
        self.assertFalse(self.theorem["strict_physical_value_generator_exported"])
        self.assertFalse(self.theorem["qw2191_discharged"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["legacy_role_transfer_exported"])
        self.assertFalse(self.theorem["toe_closure_exported"])

    def test_rg_audit_docs_and_lay_summary(self) -> None:
        audit = self.payload["rg_non_duplication_audit"]
        self.assertEqual(audit["tool"], "rg")
        self.assertIn("new_packet", audit["patterns"])
        md_text = MD.read_text(encoding="utf-8")
        self.assertIn("Plain-language progress note", md_text)
        self.assertIn("recalculates the most important saved cells", self.theorem["lay_summary"])
        self.assertIn("P2474/S1424", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2474/S1424", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
