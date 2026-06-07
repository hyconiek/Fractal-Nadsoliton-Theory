from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from p2558_s1508_legacy_to_strict_damping_power_mean_homotopy_continuum_certificate import MD, OUT, append_doc_sections, build_payload, write_markdown


class P2558DampingPowerMeanHomotopyContinuumTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload()
        OUT.write_text(json.dumps(cls.payload, indent=2, sort_keys=True), encoding="utf-8")
        write_markdown(cls.payload)
        append_doc_sections()
        cls.theorem = cls.payload["legacy_to_strict_damping_power_mean_homotopy_continuum_certificate"]["theorem_export"]

    def test_power_mean_continuum_witness(self) -> None:
        self.assertEqual(self.payload["packet_id"], "P2558")
        self.assertEqual(self.payload["stage_id"], "S1508")
        self.assertTrue(self.theorem["p2556_homotopy_nonuniqueness_inherited"])
        self.assertTrue(self.theorem["p2557_metric_dependence_inherited"])
        self.assertEqual(self.theorem["power_mean_homotopy_row_count"], 22)
        self.assertEqual(self.theorem["q_value_count"], 5)
        self.assertEqual(self.theorem["power_mean_q_values_audited"], [-2.0, -1.0, 0.0, 1.0, 2.0])
        self.assertTrue(self.theorem["all_power_mean_homotopies_share_endpoints"])
        self.assertTrue(self.theorem["all_power_mean_homotopies_share_endpoint_log_primitive"])
        self.assertTrue(self.theorem["midpoint_source_spread_positive_for_all_rows"])
        self.assertTrue(self.theorem["finite_q_sample_witnesses_continuum_nonuniqueness"])

    def test_negative_controls_and_docs(self) -> None:
        self.assertFalse(self.theorem["power_mean_homotopy_selector_source_exported"])
        self.assertFalse(self.theorem["unique_damping_homotopy_source_exported"])
        self.assertFalse(self.theorem["damping_compression_bridge_component_ready"])
        self.assertFalse(self.theorem["full_bridge_theorem_exported"])
        self.assertFalse(self.theorem["role_transfer_theorem_exported"])
        self.assertFalse(self.theorem["toe_closure_claimed"])
        self.assertIn("power-mean parameter q", self.theorem["recommended_next_honest_step"])
        self.assertTrue(all(self.payload["gatekeeper_checks"].values()))
        self.assertEqual(self.payload["rg_non_duplication_audit"]["tool"], "rg")
        self.assertIn("P2558/S1508", MD.read_text(encoding="utf-8"))
        self.assertIn("P2558/S1508", (ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md").read_text(encoding="utf-8"))
        self.assertIn("P2558/S1508", (ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
