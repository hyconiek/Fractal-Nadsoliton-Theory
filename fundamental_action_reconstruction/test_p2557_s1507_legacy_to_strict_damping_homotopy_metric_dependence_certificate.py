from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2557_s1507_legacy_to_strict_damping_homotopy_metric_dependence_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2557DampingHomotopyMetricDependenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["legacy_to_strict_damping_homotopy_metric_dependence_certificate"]["theorem_export"]

    def test_metric_dependence(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2557")
        self.assertEqual(self.payload["stage_id"], "S1507")
        self.assertTrue(self.theorem["p2556_homotopy_nonuniqueness_inherited"])
        self.assertEqual(self.theorem["row_count"], 22)
        self.assertTrue(self.theorem["denominator_velocity_metric_selects_linear_for_all_rows"])
        self.assertTrue(self.theorem["log_source_metric_selects_geometric_for_all_rows"])
        self.assertTrue(self.theorem["metrics_disagree_for_all_rows"])
        self.assertTrue(self.theorem["homotopy_metric_choice_is_additional_source_obligation"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["homotopy_metric_selector_source_exported"])
        self.assertFalse(self.theorem["unique_damping_homotopy_source_exported"])
        self.assertFalse(self.theorem["damping_compression_bridge_component_ready"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertIn("homotopy metric/source-density principle", self.theorem["recommended_next_honest_step"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2557/S1507", MD.read_text(encoding="utf-8"))
        self.assertIn("P2557/S1507", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2557/S1507", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
