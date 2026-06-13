from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "p2665_s1615_selector_lane_to_boundary_phase_sector_bridge_audit.py"
OUT = ROOT / "generated" / "p2665_s1615_selector_lane_to_boundary_phase_sector_bridge_audit.json"
MD = ROOT / "generated" / "p2665_s1615_selector_lane_to_boundary_phase_sector_bridge_audit.md"


class P2665SelectorLaneToBoundaryPhaseSectorBridgeAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        cls.payload = json.loads(OUT.read_text(encoding="utf-8"))

    def test_content_first_grep_includes_selector_and_no_name_only_mode(self) -> None:
        audit = self.payload["semantic_rg_antiduplication_audit"]
        self.assertIn("content-first", audit["mode"])
        self.assertIn("selector", audit["mode"])
        for key in (
            "selector_general_content",
            "source_topology_selector_content",
            "boundary_phase_sector_content",
            "chart_quotient_limit_content",
            "raw_theta_uniqueness_content",
            "nonclosure_guard_content",
        ):
            self.assertIn(key, audit["patterns"])
            self.assertGreater(audit["patterns"][key]["count"], 0)
        presence = self.payload["lane_presence_from_content_grep"]
        self.assertTrue(presence["selector_material_exists"])
        self.assertTrue(presence["source_topology_selector_material_exists"])
        self.assertTrue(presence["boundary_phase_sector_material_exists"])

    def test_upstream_p2664_consistency(self) -> None:
        upstream = self.payload["upstream_consistency"]
        self.assertTrue(upstream["p2664_nonexact_bit_carrier_still_available"])
        self.assertTrue(upstream["p2664_local_even_class_does_not_select"])
        self.assertTrue(upstream["p2664_declared_theta_can_select_but_not_internal"])
        self.assertTrue(upstream["p2664_boundary_phase_bit_target_not_exported"])

    def test_transfer_acceptance_matrix_blocks_selector_false_transfer(self) -> None:
        rows = {row["lane"]: row for row in self.payload["transfer_acceptance_rows"]}
        for lane in (
            "global_projective_selector_state",
            "source_topology_selector_witness",
            "declared_scope_or_basis_free_quotient_selector",
            "declared_theta_holonomy_source",
            "bridge_completed_boundary_phase_sector_selector_target",
        ):
            self.assertIn(lane, rows)
            self.assertFalse(rows[lane]["passes_bridge_acceptance_now"])
        self.assertEqual(rows["global_projective_selector_state"]["best_sectors"], [0, 1])
        self.assertTrue(rows["declared_theta_holonomy_source"]["selects_nonexact_sector_one"])
        self.assertFalse(rows["declared_theta_holonomy_source"]["strict_non_declared_boundary_phase_provenance"])

    def test_no_closure_and_docs_updated(self) -> None:
        decision = self.payload["closure_decision"]
        self.assertEqual(decision["passing_selector_to_boundary_phase_bridge_lanes"], [])
        self.assertFalse(decision["selector_lane_to_boundary_phase_bridge_exported_now"])
        self.assertFalse(decision["boundary_phase_sector_selector_exported_now"])
        self.assertFalse(decision["boundary_phase_bit_target_exported_now"])
        self.assertFalse(decision["unconditional_uv_unit_selected_now"])
        self.assertFalse(decision["beta_source_exported_now"])
        self.assertFalse(decision["qw2191_discharged_now"])
        self.assertFalse(decision["role_bearing_ltotal_now"])
        self.assertFalse(decision["toe_closure_now"])
        self.assertTrue(all(value is False for value in self.payload["negative_export_flags"].values()))
        self.assertIn("P2665/S1615", MD.read_text(encoding="utf-8"))
        self.assertIn("P2665/S1615", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2665/S1615", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
