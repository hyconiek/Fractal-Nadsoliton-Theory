from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2668_s1618_orientation_datum_to_boundary_cycle_cross_invariant_audit.py"
OUT = ROOT / "generated" / "p2668_s1618_orientation_datum_to_boundary_cycle_cross_invariant_audit.json"
MD = ROOT / "generated" / "p2668_s1618_orientation_datum_to_boundary_cycle_cross_invariant_audit.md"


class P2668OrientationDatumToBoundaryCycleCrossInvariantAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_covers_orientation_and_cross_invariant(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        for key in (
            "orientation_datum_content",
            "boundary_cycle_cross_content",
            "lane_scope_limits_content",
            "selector_pair12_content",
            "nonclosure_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)

    def test_upstream_p2667_consistency(self) -> None:
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2667_sector_one_mapping_exists"])
        self.assertTrue(upstream["p2667_reversal_unforbidden"])
        self.assertTrue(upstream["p2667_no_canonical_orientation_map"])
        self.assertTrue(upstream["p2667_no_boundary_phase_bit_target"])

    def test_existing_orientation_sources_do_not_cross_to_boundary_cycle(self) -> None:
        witness = self.payload["cross_invariant_witness"]
        self.assertTrue(witness["any_existing_orientation_export_present"])
        self.assertFalse(witness["any_source_forbids_sector_swap"])
        self.assertFalse(witness["any_source_ties_pair2_to_sector1"])
        self.assertEqual(witness["passing_cross_invariant_sources"], [])
        rows = {row["source"]: row for row in witness["rows"]}
        self.assertTrue(rows["N546_S_sel_int_pair1_orientation"]["exports_orientation"])
        self.assertTrue(rows["N500_Shannon_axis_only_orientation"]["exports_orientation"])
        for row in rows.values():
            self.assertFalse(row["is_boundary_cycle_functor"])
            self.assertFalse(row["passes_cross_invariant_acceptance_now"])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_cross_invariant_sources"], [])
        self.assertTrue(decision["existing_orientation_export_present"])
        self.assertFalse(decision["source_forbids_sector_swap"])
        self.assertFalse(decision["source_ties_pair2_to_sector1"])
        self.assertFalse(decision["boundary_phase_bit_target_exported_now"])
        self.assertFalse(decision["unconditional_uv_unit_selected_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["qw2191_discharged_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2668/S1618", MD.read_text(encoding="utf-8"))
        self.assertIn("P2668/S1618", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2668/S1618", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
