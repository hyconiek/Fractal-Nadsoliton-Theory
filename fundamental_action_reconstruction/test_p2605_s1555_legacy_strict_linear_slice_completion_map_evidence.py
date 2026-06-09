from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2605_s1555_legacy_strict_linear_slice_completion_map_evidence import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2605LegacyStrictLinearSliceCompletionMapEvidenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["legacy_strict_linear_slice_completion_map_evidence"]["theorem_export"]
        cls.audit = cls.theorem["completion_map_audit"]
        cls.remaining = cls.theorem["post_linear_slice_remaining_bridge_matrix"]

    def test_identity_and_p2604_inheritance(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2605")
        self.assertEqual(self.payload["stage_id"], "S1555")
        self.assertIn("LINEAR_SLICE_COMPLETION_MAP_EVIDENCE", self.payload["status"])
        self.assertTrue(self.theorem["p2604_strict_damping_source_inherited"])
        self.assertTrue(self.payload["gatekeeper_checks"]["p2604_source_readiness_inherited"])

    def test_linear_slice_completion_map_exactness(self) -> None:
        self.assertTrue(self.theorem["linear_slice_completion_map_evidence_exported"])
        self.assertTrue(self.audit["linear_slice_exact_on_grid"])
        self.assertLess(self.audit["max_abs_residual_on_audit_grid"], 1e-14)
        self.assertEqual(self.audit["parameter_map_evidence"]["eta_slice"], "eta <- 1")
        accepting = [row for row in self.audit["eta_candidate_rows"] if row["exact_linear_slice_accepts"]]
        self.assertEqual(len(accepting), 1)
        self.assertEqual(accepting[0]["eta_candidate"], 1.0)

    def test_remaining_bridge_matrix(self) -> None:
        self.assertEqual(self.remaining["truth_table_row_count"], 4)
        self.assertEqual(self.remaining["accepting_row_count"], 1)
        self.assertEqual(self.remaining["remaining_missing_gate_count_after_p2605"], 2)
        self.assertFalse(self.remaining["current_role_bearing_ltotal_accepts"])
        self.assertEqual(self.remaining["current_assignment_after_p2605"]["legacy_to_strict_completion_map_evidence"], True)

    def test_scope_guards_and_docs(self) -> None:
        self.assertFalse(self.theorem["full_legacy_to_strict_completion_bridge_exported"])
        self.assertFalse(self.theorem["strict_side_residual_additions_exported"])
        self.assertFalse(self.theorem["strict_damping_role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["legacy_physical_role_transfer_exported"])
        self.assertFalse(self.theorem["role_bearing_ltotal_exported"])
        self.assertFalse(self.theorem["qw2191_discharged_by_this_evidence"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2605/S1555", MD.read_text(encoding="utf-8"))
        self.assertIn("P2605/S1555", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2605/S1555", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
